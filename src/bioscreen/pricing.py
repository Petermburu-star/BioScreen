"""
BioScreen — Actuarial risk pricing framework (v2).

Corrected structure:
    risk = screening_base × context_modulator

Screening signal is dominant. Verification, history, organism, and size
modulate it — they don't add independently.
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass


@dataclass
class OrderRisk:
    customer_verified: bool
    order_size_bp: int
    organism: str
    screening_hazard_delta: float
    screening_best_toxin_sim: float
    customer_history: float


class ActuarialPricer:
    # Base rate: probability a random order is hazardous
    BASE_HAZARD_RATE = 0.001

    # Catastrophic loss estimate per biosecurity incident (USD)
    CATASTROPHIC_LOSS_USD = 1e9

    # Organism risk multipliers (1.0 = human = highest)
    ORGANISM_RISK = {
        "human": 1.0, "animal": 0.6, "plant": 0.3,
        "environmental": 0.2, "synthetic": 0.8,
    }

    def compute_screening_base(self, order: OrderRisk) -> float:
        """
        Base risk from the screening signal alone.
        - hazard_delta of 0.00 → 0
        - hazard_delta of 0.10+ → 1.0
        - toxin_similarity above 0.90 boosts the base
        """
        delta_component = min(1.0, max(0.0, order.screening_hazard_delta / 0.10))
        sim_component = max(0.0, min(1.0, (order.screening_best_toxin_sim - 0.90) / 0.10))
        return max(delta_component, sim_component)

    def compute_context_modulator(self, order: OrderRisk) -> float:
        """
        Context modulates the screening signal.
        1.0 = neutral. >1 amplifies. <1 dampens.
        """
        modulator = 1.0

        # Verification reduces effective risk
        if order.customer_verified:
            modulator *= 0.7

        # Customer history (0 = bad, 1 = excellent)
        history_bonus = 1.0 - (order.customer_history * 0.3)
        modulator *= history_bonus

        # Organism risk
        modulator *= self.ORGANISM_RISK.get(order.organism.lower(), 0.5)

        # Size: large orders amplify (log-scaled, not linear)
        size_factor = 1.0 + min(1.0, order.order_size_bp / 10000) * 0.3
        modulator *= size_factor

        return modulator

    def compute_risk_score(self, order: OrderRisk) -> float:
        base = self.compute_screening_base(order)
        mod = self.compute_context_modulator(order)
        return float(min(1.0, max(0.0, base * mod)))

    def compute_premium(self, order: OrderRisk,
                         order_cost_usd: float = 1000.0) -> dict:
        """
        Premium = base_rate × catastrophic_loss × risk_multiplier
                  + order_cost × risk_loading

        Calibrated so that:
          - LOW:      premium ≈ 0.02% of order
          - MEDIUM:   premium ≈ 0.5% of order
          - HIGH:     premium ≈ 5% of order
          - CRITICAL: premium ≈ 25% of order
          - REJECTED: rejected outright
        """
        risk = self.compute_risk_score(order)

        # Risk-loading on order cost (the main lever)
        # Piecewise calibration for tier separation
        if risk < 0.10:
            loading = 0.0010              # 0.10% — baseline biosecurity fee
        elif risk < 0.30:
            loading = 0.001 + (risk - 0.10) * 0.025   # 0.1% to 0.6%
        elif risk < 0.55:
            loading = 0.005 + (risk - 0.30) * 0.10    # 0.5% to 3%
        elif risk < 0.80:
            loading = 0.030 + (risk - 0.55) * 0.35    # 3% to 12%
        else:
            loading = 0.120 + (risk - 0.80) * 0.65    # 12% to 25%

        premium_usd = order_cost_usd * loading

        # Tier classification
        if risk < 0.10:
            tier = "LOW"
            action = "Standard screening; base price"
        elif risk < 0.30:
            tier = "MEDIUM"
            action = "Enhanced screening required"
        elif risk < 0.55:
            tier = "HIGH"
            action = "Full verification + significant surcharge"
        elif risk < 0.80:
            tier = "CRITICAL"
            action = "Senior review required; prohibitive surcharge"
        else:
            tier = "REJECTED"
            action = "De facto rejection; institutional exemption required"

        return {
            "risk_score": round(risk, 4),
            "premium_usd": round(premium_usd, 2),
            "premium_pct_of_order": round(premium_usd / max(order_cost_usd, 1) * 100, 3),
            "tier": tier,
            "action": action,
        }

    def simulate_market(self, n_orders: int = 5000) -> pd.DataFrame:
        np.random.seed(42)
        rows = []

        for _ in range(n_orders):
            # 1% hazardous orders
            if np.random.random() < 0.01:
                delta = np.random.uniform(0.05, 0.20)
                toxin_sim = np.random.uniform(0.95, 1.00)
            else:
                delta = np.random.beta(1.5, 30) * 0.04
                toxin_sim = np.random.beta(5, 3) * 0.25 + 0.65

            order = OrderRisk(
                customer_verified=np.random.random() < 0.75,
                order_size_bp=int(np.random.exponential(1500)) + 100,
                organism=np.random.choice(
                    ["human", "animal", "plant", "environmental", "synthetic"],
                    p=[0.4, 0.2, 0.15, 0.15, 0.10]
                ),
                screening_hazard_delta=delta,
                screening_best_toxin_sim=toxin_sim,
                customer_history=np.random.beta(5, 2),
            )

            order_cost = np.random.exponential(2000) + 500
            result = self.compute_premium(order, order_cost_usd=order_cost)

            rows.append({
                **order.__dict__,
                "order_cost_usd": round(order_cost, 2),
                **result,
            })

        return pd.DataFrame(rows)
