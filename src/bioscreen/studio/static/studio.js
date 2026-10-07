const $ = (id) => document.getElementById(id);
let lastResult = null;

/* ---------- status ---------- */
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

/* ---------- sequence input ---------- */
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

/* ---------- color helpers ---------- */
function lerpColor(a, b, t) {
  return [
    Math.round(a[0] + (b[0] - a[0]) * t),
    Math.round(a[1] + (b[1] - a[1]) * t),
    Math.round(a[2] + (b[2] - a[2]) * t),
  ];
}

/* ---------- ribbon renderer with adaptive normalization ---------- */
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

  // ---- adaptive normalization: fit to this profile's own range ----
  let vmin = Infinity, vmax = -Infinity;
  for (const v of values) {
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }
  const span = Math.max(vmax - vmin, 1e-6);
  // pad the range slightly so darkest/brightest pixels are visible
  const lo = vmin - span * 0.05;
  const hi = vmax + span * 0.05;
  const range = hi - lo;

  ctx.fillStyle = "#070b14";
  ctx.fillRect(0, 0, W, H);

  for (let x = 0; x < W; x++) {
    const i0 = Math.floor((x / W) * L);
    const i1 = Math.max(i0 + 1, Math.floor(((x + 1) / W) * L));
    let v = 0;
    for (let i = i0; i < i1 && i < L; i++) {
      if (values[i] > v) v = values[i];
    }
    let t = (v - lo) / range;
    t = Math.max(0, Math.min(1, t));

    let color;
    if (mode === "toxin") {
      // dark blue -> amber (bright)
      color = lerpColor([16, 24, 40], [245, 158, 11], t);
    } else {
      // dark blue -> cyan
      color = lerpColor([16, 24, 40], [34, 211, 238], t);
    }
    ctx.fillStyle = "rgb(" + color[0] + "," + color[1] + "," + color[2] + ")";
    ctx.fillRect(x, 0, 1, H);
  }

  // ---- residue ruler: ticks every 100 residues ----
  ctx.strokeStyle = "rgba(148,163,184,.20)";
  ctx.lineWidth = 1;
  ctx.fillStyle = "rgba(148,163,184,.55)";
  ctx.font = "9px monospace";
  for (let i = 100; i < L; i += 100) {
    const x = Math.round((i / L) * W);
    ctx.beginPath();
    ctx.moveTo(x, H - 6);
    ctx.lineTo(x, H);
    ctx.stroke();
    if (i % 500 === 0) {
      ctx.fillText(String(i), x + 2, H - 8);
    }
  }
}

/* ---------- divergence renderer with mean baseline ---------- */
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

  // compute per-residue divergence and its mean
  const diffs = new Array(L);
  let sum = 0, absMax = 0;
  for (let i = 0; i < L; i++) {
    const d = toxin[i] - safeValues[i];
    diffs[i] = d;
    sum += d;
    if (Math.abs(d) > absMax) absMax = Math.abs(d);
  }
  const mean = sum / L;
  // scale so the largest deviation from the mean maps to the canvas half-height
  const scale = Math.max(absMax, 0.005);

  ctx.fillStyle = "#070b14";
  ctx.fillRect(0, 0, W, H);

  // ---- centre reference line (zero) ----
  ctx.strokeStyle = "rgba(148,163,184,.15)";
  ctx.beginPath();
  ctx.moveTo(0, mid);
  ctx.lineTo(W, mid);
  ctx.stroke();

  // ---- mean baseline ----
  const meanY = mid - (mean / scale) * (mid - 4);
  ctx.strokeStyle = "rgba(148,163,184,.55)";
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(0, meanY);
  ctx.lineTo(W, meanY);
  ctx.stroke();
  ctx.setLineDash([]);

  // ---- columns ----
  for (let x = 0; x < W; x++) {
    const i0 = Math.floor((x / W) * L);
    const i1 = Math.max(i0 + 1, Math.floor(((x + 1) / W) * L));
    // pick the value with largest |deviation from mean| in this pixel column
    let bestDev = 0, bestD = 0;
    for (let i = i0; i < i1 && i < L; i++) {
      const dev = diffs[i] - mean;
      if (Math.abs(dev) > Math.abs(bestDev)) {
        bestDev = dev;
        bestD = diffs[i];
      }
    }
    const y = mid - (bestDev / scale) * (mid - 4);
    const top = Math.min(y, meanY);
    const bot = Math.max(y, meanY);
    if (bestD - mean >= 0) {
      ctx.fillStyle = "rgba(245,158,11,.85)";
    } else {
      ctx.fillStyle = "rgba(34,211,238,.85)";
    }
    ctx.fillRect(x, top, 1, Math.max(1, bot - top));
  }
}

/* ---------- gate bar updater ---------- */
function updateGate(containerId, value, threshold, format) {
  const el = $(containerId);
  if (!el) return;
  const passed = value >= threshold;
  el.classList.toggle("pass", passed);
  el.classList.toggle("fail", !passed);

  const fill = el.querySelector(".gate-fill");
  const thr = el.querySelector(".gate-threshold");
  const val = el.querySelector(".gate-val");

  // Scale the bar so [0 .. max(threshold*2, value*1.1)] fills the width
  const max = Math.max(threshold * 2, value * 1.1, 0.01);
  fill.style.width = Math.min(100, (value / max) * 100) + "%";
  thr.style.left = Math.min(100, (threshold / max) * 100) + "%";
  val.textContent = format(value);
}

/* ---------- screen ---------- */
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
  $("fingerprint-panel").classList.remove("hidden");
  $("verdict-panel").classList.remove("hidden");
  $("pricing-panel").classList.remove("hidden");
  renderFingerprint(j);

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

  // gate bars
  updateGate("gate-delta", j.hazard_delta, j.threshold_delta,
             v => v.toFixed(4));
  updateGate("gate-floor", j.best_toxin_score, j.threshold_floor,
             v => v.toFixed(4));
}

/* ---------- pricing ---------- */
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

  // provenance line
  const orderCost = parseFloat($("in-cost").value) || 1000;
  const prov = $("out-provenance");
  if (prov) {
    prov.textContent = "$" + orderCost.toFixed(0)
      + "  x  " + j.premium_pct_of_order.toFixed(3) + "%"
      + "  =  $" + j.premium_usd.toFixed(2);
  }
}

/* ---------- bind controls ---------- */
["in-verified", "in-organism", "in-size", "in-cost"].forEach((id) => {
  $(id).addEventListener("input", renderPricing);
});
$("in-history").addEventListener("input", () => {
  $("in-history-val").textContent = parseFloat($("in-history").value).toFixed(2);
  renderPricing();
});

/* ---------- copy json ---------- */
const copyBtn = $("btn-copy");
if (copyBtn) {
  copyBtn.addEventListener("click", async () => {
    if (!lastResult) return;
    const text = JSON.stringify(lastResult, null, 2);
    try {
      await navigator.clipboard.writeText(text);
      const orig = copyBtn.textContent;
      copyBtn.textContent = "Copied";
      setTimeout(() => (copyBtn.textContent = orig), 1400);
    } catch (e) {
      alert("Clipboard unavailable");
    }
  });
}

/* ---------- resize ---------- */
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (lastResult) {
      renderRibbon("toxin-ribbon", lastResult.toxin_profile, "toxin");
      renderRibbon("safe-ribbon", lastResult.safe_profile, "safe");
      renderDivergence("diverge-ribbon",
                        lastResult.toxin_profile,
                        lastResult.safe_profile);
    }
  }, 150);
});

warmup();


/* ---------- fingerprint panel rendering ---------- */
function renderFingerprint(j) {
  const panel = $("fingerprint-panel");
  if (!panel) return;
  panel.classList.remove("hidden");

  const hero = $("fp-hero");
  const flagged = j.fingerprint_flagged;
  hero.className = "fp-hero " + (flagged ? "" : "clear");

  $("fp-glyph").textContent = flagged ? "◆" : "◇";
  $("fp-name").textContent = j.fingerprint_best + "  ·  " + j.fingerprint_best_accession;
  $("fp-score").textContent = j.fingerprint_score.toFixed(4);
  $("fp-top-k").textContent = j.fingerprint_top_k.toFixed(4);
  $("fp-order").textContent = j.fingerprint_order.toFixed(3);

  // Ranking table
  const rank = $("fp-ranking");
  rank.innerHTML = "";
  j.fingerprint_ranking.forEach((row, idx) => {
    const el = document.createElement("div");
    el.className = "fp-rank-row";
    el.innerHTML =
      '<span class="nm">' + (idx + 1) + '. ' + row.name + '</span>' +
      '<div class="fp-rank-bar"><span style="width:' + Math.min(100, row.score * 100) + '%"></span></div>' +
      '<span class="sc">' + row.score.toFixed(4) + '</span>';
    rank.appendChild(el);
  });
}


/* ---------- EvoDiff variant chips ---------- */
let _evodiffCache = null;

async function _loadEvoDiffVariants() {
  if (_evodiffCache) return _evodiffCache;
  const r = await fetch("/api/evodiff_variants");
  if (!r.ok) throw new Error("evodiff_variants fetch failed");
  _evodiffCache = await r.json();
  return _evodiffCache;
}

document.querySelectorAll(".chip.evo").forEach((btn) => {
  btn.addEventListener("click", async () => {
    const idx = parseInt(btn.dataset.evo, 10);
    const original = btn.textContent;
    btn.textContent = "...";
    try {
      const data = await _loadEvoDiffVariants();
      const v = data.variants[idx];
      if (!v) {
        alert("Variant " + idx + " not found");
        return;
      }
      $("seq-input").value = v.sequence;
      if (typeof updateSeqMeta === "function") updateSeqMeta();
      // Show the pre-computed scores in a small banner
      const note = document.createElement("div");
      note.className = "evo-banner";
      note.innerHTML =
        '<strong>EvoDiff variant v' + idx + '</strong>  ' +
        'identity to native ricin = ' + (v.sequence_identity_to_native * 100).toFixed(2) + '%  ·  ' +
        'BLAST score = ' + v.blast_score.toFixed(3) + ' (evaded)  ·  ' +
        'fingerprint = ' + v.fingerprint_score.toFixed(3) + ' (caught)  ·  ' +
        '<em>Click Screen to re-validate live</em>';
      const existing = document.querySelector(".evo-banner");
      if (existing) existing.remove();
      const wrap = $("seq-input").parentElement;
      wrap.appendChild(note);
    } catch (e) {
      alert("Could not load variant: " + e.message);
    }
    btn.textContent = original;
  });
});
