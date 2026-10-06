const $ = (id) => document.getElementById(id);
let lastResult = null;

async function pollStatus() {
  try {
    const r = await fetch("/api/status");
    const j = await r.json();
    if (j.status === "ready") {
      $("status-pip").className = "pip ready";
      $("status-text").textContent = "engine ready - " + j.toxins + " toxins - " + j.safes + " safes";
      $("btn-screen").disabled = false;
    } else if (j.status === "loading" || j.status === "idle") {
      $("status-pip").className = "pip loading";
      $("status-text").textContent = "loading engine (ESM-C 600M)...";
      setTimeout(pollStatus, 1500);
    } else {
      $("status-pip").className = "pip error";
      $("status-text").textContent = "error: " + (j.error || "unknown");
    }
  } catch (e) {
    $("status-pip").className = "pip error";
    $("status-text").textContent = "backend offline";
  }
}

async function warmup() {
  await fetch("/api/warmup", { method: "POST" });
  pollStatus();
}

function updateSeqMeta() {
  const s = $("seq-input").value.replace(/\s/g, "").toUpperCase();
  $("seq-len").textContent = s.length + " aa";
  if (s.length > 0) {
    const counts = {};
    for (const aa of s) counts[aa] = (counts[aa] || 0) + 1;
    const top = Object.entries(counts).sort((a, b) => b[1] - a[1]).slice(0, 4);
    $("seq-comp").textContent = top
      .map(([k, v]) => k + " " + Math.round((v / s.length) * 100) + "%")
      .join(" - ");
  } else {
    $("seq-comp").textContent = "-";
  }
}
$("seq-input").addEventListener("input", updateSeqMeta);

document.querySelectorAll(".chip").forEach((btn) => {
  btn.addEventListener("click", async () => {
    const name = btn.dataset.ex;
    const original = btn.textContent;
    btn.textContent = "...";
    try {
      const r = await fetch("/api/example/" + name);
      const j = await r.json();
      if (j.sequence) {
        $("seq-input").value = j.sequence;
        updateSeqMeta();
      } else {
        alert(j.error || "Could not load example");
      }
    } catch (e) {
      alert("UniProt fetch failed");
    }
    btn.textContent = original;
  });
});

function lerpColor(a, b, t) {
  return [
    Math.round(a[0] + (b[0] - a[0]) * t),
    Math.round(a[1] + (b[1] - a[1]) * t),
    Math.round(a[2] + (b[2] - a[2]) * t),
  ];
}

function renderRibbon(canvasId, values, mode) {
  const c = $(canvasId);
  if (!c) return;
  const dpr = window.devicePixelRatio || 1;
  const W = c.clientWidth;
  const H = c.clientHeight;
  c.width = W * dpr;
  c.height = H * dpr;
  const ctx = c.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, W, H);

  const L = values.length;
  if (L === 0) return;

  ctx.fillStyle = "#070b14";
  ctx.fillRect(0, 0, W, H);

  for (let x = 0; x < W; x++) {
    const i0 = Math.floor((x / W) * L);
    const i1 = Math.max(i0 + 1, Math.floor(((x + 1) / W) * L));
    let v = 0;
    for (let i = i0; i < i1 && i < L; i++) {
      if (values[i] > v) v = values[i];
    }
    let t = (v - 0.60) / 0.40;
    t = Math.max(0, Math.min(1, t));

    let color;
    if (mode === "toxin") color = lerpColor([13, 21, 38], [245, 158, 11], t);
    else color = lerpColor([13, 21, 38], [34, 211, 238], t);

    ctx.fillStyle = "rgb(" + color[0] + "," + color[1] + "," + color[2] + ")";
    ctx.fillRect(x, 0, 1, H);
  }

  ctx.strokeStyle = "rgba(255,255,255,.05)";
  ctx.lineWidth = 1;
  for (let i = 100; i < L; i += 100) {
    const x = Math.round((i / L) * W);
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, H);
    ctx.stroke();
  }
}

function renderDivergence(canvasId, toxin, safeValues) {
  const c = $(canvasId);
  if (!c) return;
  const dpr = window.devicePixelRatio || 1;
  const W = c.clientWidth;
  const H = c.clientHeight;
  c.width = W * dpr;
  c.height = H * dpr;
  const ctx = c.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, W, H);

  const L = toxin.length;
  if (L === 0) return;
  const mid = H / 2;

  ctx.fillStyle = "#070b14";
  ctx.fillRect(0, 0, W, H);

  ctx.strokeStyle = "rgba(148,163,184,.25)";
  ctx.beginPath();
  ctx.moveTo(0, mid);
  ctx.lineTo(W, mid);
  ctx.stroke();

  for (let x = 0; x < W; x++) {
    const i0 = Math.floor((x / W) * L);
    const i1 = Math.max(i0 + 1, Math.floor(((x + 1) / W) * L));
    let d = 0;
    for (let i = i0; i < i1 && i < L; i++) {
      const v = toxin[i] - safeValues[i];
      if (Math.abs(v) > Math.abs(d)) d = v;
    }
    const bh = Math.min(mid, Math.abs((d / 0.10) * mid));
    if (d >= 0) {
      ctx.fillStyle = "rgb(245,158,11)";
      ctx.fillRect(x, mid - bh, 1, bh);
    } else {
      ctx.fillStyle = "rgb(34,211,238)";
      ctx.fillRect(x, mid, 1, bh);
    }
  }
}

$("btn-screen").addEventListener("click", async () => {
  const seq = $("seq-input").value.replace(/\s/g, "").toUpperCase();
  if (seq.length < 20) return;

  const btn = $("btn-screen");
  btn.classList.add("busy");
  btn.disabled = true;
  btn.innerHTML = '<span class="btn-icon">&#9680;</span> screening...';

  try {
    const r = await fetch("/api/screen", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence: seq }),
    });
    const j = await r.json();
    if (j.error) throw new Error(j.error);
    lastResult = j;
    renderResult(j);
    await renderPricing();
  } catch (e) {
    alert("Screen failed: " + e.message);
  } finally {
    btn.classList.remove("busy");
    btn.disabled = false;
    btn.innerHTML = '<span class="btn-icon">&#9656;</span> Screen';
  }
});

function renderResult(j) {
  $("instrument").classList.remove("hidden");
  $("verdict-panel").classList.remove("hidden");
  $("pricing-panel").classList.remove("hidden");

  $("toxin-ref").textContent = j.best_toxin_name + " - " + j.best_toxin_acc;
  $("toxin-score").textContent = j.best_toxin_score.toFixed(4);
  $("safe-ref").textContent = j.best_safe_name + " - " + j.best_safe_acc;
  $("safe-score").textContent = j.best_safe_score.toFixed(4);

  renderRibbon("toxin-ribbon", j.toxin_profile, "toxin");
  renderRibbon("safe-ribbon", j.safe_profile, "safe");
  renderDivergence("diverge-ribbon", j.toxin_profile, j.safe_profile);

  const hero = $("verdict-hero");
  hero.className = "verdict-hero " + (j.flagged ? "flagged" : "clear");
  $("verdict-glyph").textContent = j.flagged ? "\u26A0" : "\u2713";
  $("verdict-word").textContent = j.flagged ? "Flagged" : "Clear";
  $("verdict-line").textContent = j.flagged
    ? "This sequence carries a functional signature of a controlled toxin."
    : "No functional toxin signature detected.";

  $("thr-delta").textContent = j.threshold_delta.toFixed(2);
  $("thr-floor").textContent = j.threshold_floor.toFixed(2);
  $("gate-delta-val").textContent = j.hazard_delta.toFixed(4);
  $("gate-floor-val").textContent = j.best_toxin_score.toFixed(4);
  $("gate-delta").classList.toggle("open", j.gate_delta);
  $("gate-floor").classList.toggle("open", j.gate_floor);
}

async function renderPricing() {
  if (!lastResult) return;
  const body = {
    verified: $("in-verified").checked,
    history: parseFloat($("in-history").value),
    organism: $("in-organism").value,
    size_bp: parseInt($("in-size").value) || 2000,
    order_cost: parseFloat($("in-cost").value) || 1000,
    hazard_delta: lastResult.hazard_delta,
    toxin_sim: lastResult.best_toxin_score,
  };
  const r = await fetch("/api/price", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const j = await r.json();

  const ro = document.querySelector(".readout");
  ro.className = "readout tier-" + j.tier;

  $("out-tier").textContent = j.tier;
  $("out-premium").textContent = "$" + j.premium_usd.toFixed(2);
  $("out-risk").textContent = j.risk_score.toFixed(4);
  $("out-loading").textContent = j.premium_pct_of_order.toFixed(3) + "%";
  $("out-action").textContent = j.action;
}

["in-verified", "in-organism", "in-size", "in-cost"].forEach((id) => {
  $(id).addEventListener("input", renderPricing);
});
$("in-history").addEventListener("input", () => {
  $("in-history-val").textContent = parseFloat($("in-history").value).toFixed(2);
  renderPricing();
});

let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (lastResult) {
      renderRibbon("toxin-ribbon", lastResult.toxin_profile, "toxin");
      renderRibbon("safe-ribbon", lastResult.safe_profile, "safe");
      renderDivergence("diverge-ribbon", lastResult.toxin_profile, lastResult.safe_profile);
    }
  }, 150);
});

warmup();
