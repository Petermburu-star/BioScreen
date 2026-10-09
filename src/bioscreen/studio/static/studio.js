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
    } else if (j.status === "loading") {
      $("status-pip").className = "pip loading";
      const p = j.progress || {};
      const secs = p.elapsed_s !== undefined ? ` [${p.elapsed_s}s]` : "";
      const stage = p.stage ? p.stage.toUpperCase() : "loading";
      const detail = p.detail || "loading engine";
      $("status-text").textContent = stage + secs + " · " + detail;
      setTimeout(pollStatus, 1000);
    } else if (j.status === "idle") {
      $("status-pip").className = "pip loading";
      $("status-text").textContent = "waiting for engine to start...";
      setTimeout(pollStatus, 1000);
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
  renderTrace(j, "screen");

  $("contrastive-panel").classList.remove("hidden");
  $("instrument").classList.remove("hidden");
  $("fingerprint-panel").classList.remove("hidden");
  $("verdict-panel").classList.remove("hidden");
  $("pricing-panel").classList.remove("hidden");
  renderContrastive(j);
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

  // Add verdict reason line
  const heroCopy = hero.querySelector(".verdict-copy");
  if (heroCopy) {
    let reasonEl = heroCopy.querySelector(".verdict-reason");
    if (!reasonEl) {
      reasonEl = document.createElement("div");
      reasonEl.className = "verdict-reason";
      heroCopy.appendChild(reasonEl);
    }
    reasonEl.textContent = "Signal: " + (j.verdict_reason || "—");
  }

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


/* ============================================================
   TAB SWITCHING
   ============================================================ */
document.querySelectorAll(".tab").forEach((btn) => {
  btn.addEventListener("click", () => {
    const target = btn.dataset.tab;
    document.querySelectorAll(".tab").forEach((b) => b.classList.toggle("active", b === btn));
    document.querySelectorAll(".tab-panel").forEach((p) => {
      p.classList.toggle("active", p.dataset.panel === target);
    });
  });
});


/* ============================================================
   TRACE RENDERER
   ============================================================ */
let _lastResult = null;

function renderTrace(result) {
  const container = $("trace-content");
  if (!container) return;
  container.innerHTML = "";

  if (!result || !result.trace || result.trace.length === 0) {
    container.innerHTML = '<p class="trace-empty">No trace available for this request.</p>';
    return;
  }

  // Determine whether this trace has variant tags (Generate tab) or not (Screen tab)
  const hasVariants = result.trace.some((s) => s.variant !== undefined && s.variant !== null);

  // Header — compute total defensively
  const totalMs = (result.total_elapsed_ms !== undefined && result.total_elapsed_ms !== null)
    ? result.total_elapsed_ms
    : result.trace.reduce((sum, s) => sum + (s.elapsed_ms || 0), 0);

  const hdr = document.createElement("div");
  hdr.className = "trace-header";
  hdr.innerHTML =
    '<span class="th-left">' + result.trace.length + ' steps executed</span>' +
    '<span class="th-total">total: ' + totalMs.toFixed(1) + ' ms</span>';
  container.appendChild(hdr);

  // Steps
  result.trace.forEach((step) => {
    const el = document.createElement("div");
    el.className = "trace-step";
    el.dataset.step = step.step;

    let details = "";
    if (step.details && Object.keys(step.details).length > 0) {
      const lines = Object.entries(step.details)
        .map(([k, v]) => "  " + k + ": " + (typeof v === "object" ? JSON.stringify(v) : v))
        .join("\n");
      details = '<div class="ts-details"><pre>' + lines + '</pre></div>';
    }

    const variantBadge = (hasVariants && step.variant !== null && step.variant !== undefined)
      ? '<span class="ts-variant">v' + step.variant + '</span>'
      : '';

    el.innerHTML =
      '<div class="ts-phase">' + (step.phase || "step") + ' ' + variantBadge + '</div>' +
      '<div class="ts-name">' + step.name + '</div>' +
      '<div class="ts-time">' + step.elapsed_ms.toFixed(2) + ' ms</div>' +
      '<div class="ts-io">' +
        '<span class="io-k">input</span><span class="io-v">' + step.input + '</span>' +
        '<span class="io-k">output</span><span class="io-v">' + step.output + '</span>' +
      '</div>' +
      details;

    container.appendChild(el);
  });
}


/* ============================================================
   SOURCE VIEWER (Reproduce tab)
   ============================================================ */
let _sourceCache = null;

async function loadSources() {
  if (_sourceCache) return _sourceCache;
  const r = await fetch("/api/source");
  _sourceCache = await r.json();
  return _sourceCache;
}

document.querySelectorAll(".rp-src-btn").forEach((btn) => {
  btn.addEventListener("click", async () => {
    const file = btn.dataset.file;
    document.querySelectorAll(".rp-src-btn").forEach((b) => b.classList.toggle("active", b === btn));
    const sources = await loadSources();
    const view = $("source-view");
    view.textContent = sources[file] || "// file not found";
  });
});


/* ============================================================
   HOOK INTO SCREEN RESULT
   ============================================================ */
// Wrap the existing renderResult to also feed the trace panel
if (typeof window._origRenderResult === "undefined") {
  window._origRenderResult = window.renderResult;
}

/* ============================================================
   ADVERSARIAL TAB
   ============================================================ */
let _advLoaded = false;

async function loadAdversarial() {
  if (_advLoaded) return;
  const tbody = $("adv-tbody");
  if (!tbody) return;

  try {
    const r = await fetch("/api/adversarial");
    if (!r.ok) {
      tbody.innerHTML = '<tr><td colspan="7" class="adv-loading">' +
        'Results file not found. Run the 16-toxin validation first.</td></tr>';
      return;
    }
    const data = await r.json();
    renderAdversarial(data);
    _advLoaded = true;
  } catch (e) {
    tbody.innerHTML = '<tr><td colspan="7" class="adv-loading">' +
      'Load failed: ' + e.message + '</td></tr>';
  }
}

function renderAdversarial(data) {
  $("adv-n-toxins").textContent = data.n_toxins;
  $("adv-n-variants").textContent = data.n_variants;

  const blPct = ((data.blast_evaded / data.n_variants) * 100).toFixed(0);
  $("adv-blast-evaded").textContent = data.blast_evaded + " / " + data.n_variants + " (" + blPct + "%)";

  const fpPct = ((data.fingerprint_caught / data.n_variants) * 100).toFixed(0);
  $("adv-fp-caught").textContent = data.fingerprint_caught + " / " + data.n_variants + " (" + fpPct + "%)";

  const tbody = $("adv-tbody");
  tbody.innerHTML = "";

  data.toxins.forEach((tox) => {
    tox.variants.forEach((v, idx) => {
      const tr = document.createElement("tr");
      const isFirst = idx === 0;
      const blastCls = v.blast_verdict === "EVADED" ? "evaded" : "blast-caught";
      const fpCls = v.fingerprint_verdict === "CAUGHT" ? "caught" : "clear";

      tr.innerHTML =
        '<td class="adv-toxin-name">' + (isFirst ? tox.name : '') + '</td>' +
        '<td class="adv-variant-idx">v' + v.idx + '</td>' +
        '<td>' + (v.identity * 100).toFixed(2) + '%</td>' +
        '<td>' + v.blast.toFixed(3) + '</td>' +
        '<td>' + v.fingerprint.toFixed(3) + '</td>' +
        '<td><span class="adv-verdict ' + blastCls + '">' + v.blast_verdict + '</span></td>' +
        '<td><span class="adv-verdict ' + fpCls + '">' + v.fingerprint_verdict + '</span></td>';

      tbody.appendChild(tr);
    });
  });
}

// Hook: load when tab is clicked
document.querySelectorAll(".tab").forEach((btn) => {
  if (btn.dataset.tab === "adversarial") {
    btn.addEventListener("click", () => setTimeout(loadAdversarial, 50));
  }
});

// Also preload silently after a short delay so first click is instant
setTimeout(() => { if (!_advLoaded) loadAdversarial(); }, 1500);


/* ============================================================
   GENERATE TAB
   ============================================================ */
let _toxinsLoaded = false;
let _pollTimer = null;

async function loadToxins() {
  if (_toxinsLoaded) return;
  try {
    const r = await fetch("/api/toxins");
    const toxins = await r.json();
    const sel = $("gen-toxin");
    if (!sel) return;
    sel.innerHTML = "";
    toxins.forEach((t) => {
      const opt = document.createElement("option");
      opt.value = t.accession;
      opt.textContent = t.name + "  (" + t.mechanism + ")";
      sel.appendChild(opt);
    });
    _toxinsLoaded = true;
  } catch (e) {
    console.error("Failed to load toxins:", e);
  }
}

document.addEventListener("DOMContentLoaded", () => {
  const startBtn = $("gen-start");
  if (!startBtn) return;

  startBtn.addEventListener("click", async () => {
    const acc = $("gen-toxin").value;
    const mode = $("gen-mode").value;
    const n = parseInt($("gen-n").value) || 3;
    if (!acc) return;

    const btn = $("gen-start");
    btn.disabled = true;
    btn.textContent = "Generating...";

    const status = $("gen-status");
    status.style.display = "block";
    $("gen-status-stage").textContent = "queued";
    $("gen-status-pct").textContent = "0%";
    $("gen-bar-fill").style.width = "0%";
    $("gen-note").textContent = mode === "aggressive"
      ? "Aggressive mode redesigns the entire region. On CPU: 1–3 minutes per variant."
      : "Preserved mode protects the catalytic window. On CPU: 1–3 minutes per variant.";
    $("gen-results").innerHTML = "";

    try {
      const r = await fetch("/api/generate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ accession: acc, mode: mode, n_variants: n }),
      });
      const j = await r.json();
      if (j.error) throw new Error(j.error);
      pollJob(j.job_id);
    } catch (e) {
      $("gen-status-stage").textContent = "error: " + e.message;
      btn.disabled = false;
      btn.textContent = "Generate live";
    }
  });
});

function pollJob(jobId) {
  if (_pollTimer) clearInterval(_pollTimer);
  _pollTimer = setInterval(async () => {
    try {
      const r = await fetch("/api/generate/" + jobId);
      const job = await r.json();

      if (job.status === "ready") {
        clearInterval(_pollTimer);
        $("gen-status-stage").textContent = "done in " + (job.elapsed_s || "?") + " s";
        $("gen-status-pct").textContent = "100%";
        $("gen-bar-fill").style.width = "100%";
        $("gen-start").disabled = false;
        $("gen-start").textContent = "Generate live";
        renderGenerateResult(job.result);
      } else if (job.status === "error") {
        clearInterval(_pollTimer);
        $("gen-status-stage").textContent = "error: " + (job.error || "unknown");
        $("gen-start").disabled = false;
        $("gen-start").textContent = "Generate live";
      } else {
        $("gen-status-stage").textContent = job.stage || job.status;
        const pct = Math.round((job.progress || 0) * 100);
        $("gen-status-pct").textContent = pct + "%";
        $("gen-bar-fill").style.width = pct + "%";
      }
    } catch (e) {
      console.error("Poll error:", e);
    }
  }, 2000);
}



/* ============================================================
   PIPELINE TRACE — generate-tab routing
   ============================================================ */
function showGenerateTrace(result) {
  // result is the full job.result containing .trace and .variants
  if (!result || !result.trace) return;
  renderTrace(result);
  // Switch to the Pipeline trace tab automatically
  document.querySelectorAll(".tab").forEach((b) => {
    b.classList.toggle("active", b.dataset.tab === "trace");
  });
  document.querySelectorAll(".tab-panel").forEach((p) => {
    p.classList.toggle("active", p.dataset.panel === "trace");
  });
}

function renderGenerateResult(result) {
  const el = $("gen-results");
  if (!el) return;
  el.innerHTML = "";

  // Header
  const head = document.createElement("div");
  head.className = "gen-cross";
  head.innerHTML =
    '<div class="gen-cross-head">' +
      result.name + ' &middot; ' + result.mode + ' mode &middot; ' +
      result.region_length + ' aa region &middot; ' +
      result.n_catalytic + ' catalytic residues' +
    '</div>';
  el.appendChild(head);

  // Per-variant table
  const tbl = document.createElement("table");
  tbl.className = "gen-result-table";
  tbl.id = "gen-result-tbl";
  tbl.innerHTML =
    '<thead><tr>' +
      '<th>Variant</th>' +
      '<th>Identity</th>' +
      '<th>BLAST</th>' +
      '<th>Contrastive</th>' +
      '<th>Fingerprint</th>' +
      '<th>Verdict</th>' +
      '<th>Signal</th>' +
      '<th>Tier / Premium</th>' +
      '<th>Actions</th>' +
    '</tr></thead><tbody></tbody>';

  const tbody = tbl.querySelector("tbody");
  result.variants.forEach((v) => {
    const tr = document.createElement("tr");
    tr.className = "gen-variant-row";

    const blCls = v.blast_evaded ? "evaded" : "blast-caught";
    const fpCls = v.fingerprint_caught ? "caught" : "clear";
    const ctrFlagged = v.contrastive_flagged || false;
    const ctrCls = ctrFlagged ? "caught" : "clear";
    const ctrLabel = ctrFlagged ? "FLAG" : "clear";

    const verdict = v.verdict || "—";
    const verdictCls = verdict === "FLAGGED" ? "evaded" : "caught";
    const tier = (v.pricing && v.pricing.tier) ? v.pricing.tier : "—";
    const premium = (v.pricing && v.pricing.premium_usd !== undefined)
      ? "$" + v.pricing.premium_usd.toFixed(2) : "—";

    const ctrDelta = (v.contrastive_delta !== undefined && v.contrastive_delta !== null)
      ? (v.contrastive_delta >= 0 ? "+" : "") + v.contrastive_delta.toFixed(4)
      : "—";

    const reason = v.verdict_reason || "—";
    let sigClass = "clear", sigLabel = "—";
    if (reason.indexOf("Both") === 0) { sigClass = "caught"; sigLabel = "BOTH"; }
    else if (reason.indexOf("Fingerprint") === 0) { sigClass = "caught"; sigLabel = "FP"; }
    else if (reason.indexOf("Contrastive") === 0) { sigClass = "caught"; sigLabel = "CTR"; }
    else if (reason.indexOf("Neither") === 0) { sigClass = "clear"; sigLabel = "NONE"; }

    // Combined cell: score + verdict for each method
    const blastCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + v.blast.toFixed(3) + '</span>' +
        '<span class="gen-badge ' + blCls + '">' + (v.blast_evaded ? 'EVADED' : 'caught') + '</span>' +
      '</div>';
    const ctrCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + ctrDelta + '</span>' +
        '<span class="gen-badge ' + ctrCls + '">' + ctrLabel + '</span>' +
      '</div>';
    const fpCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + (v.fingerprint !== null ? v.fingerprint.toFixed(3) : '—') + '</span>' +
        '<span class="gen-badge ' + fpCls + '">' + (v.fingerprint_caught ? 'CAUGHT' : 'clear') + '</span>' +
      '</div>';
    const tierCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + tier + '</span>' +
        '<span class="gen-num-muted">' + premium + '</span>' +
      '</div>';

    tr.innerHTML =
      '<td>v' + v.idx + '</td>' +
      '<td>' + (v.identity * 100).toFixed(2) + '%</td>' +
      '<td>' + blastCell + '</td>' +
      '<td>' + ctrCell + '</td>' +
      '<td>' + fpCell + '</td>' +
      '<td><span class="gen-badge ' + verdictCls + '">' + verdict + '</span></td>' +
      '<td><span class="gen-badge ' + sigClass + '">' + sigLabel + '</span></td>' +
      '<td>' + tierCell + '</td>' +
      '<td>' +
        '<button class="gen-mini" data-vidx="' + v.idx + '" data-action="toggle-seq">seq</button> ' +
        '<button class="gen-mini primary" data-vidx="' + v.idx + '" data-action="send-to-screen">→ Screen</button>' +
      '</td>';
    tbody.appendChild(tr);

    // Hidden sequence row
    const seqTr = document.createElement("tr");
    seqTr.className = "gen-seq-row";
    seqTr.dataset.seqFor = v.idx;
    seqTr.style.display = "none";
    seqTr.innerHTML =
      '<td colspan="9">' +
        '<pre class="gen-seq-block">' + (v.sequence || "no sequence stored") + '</pre>' +
      '</td>';
    tbody.appendChild(seqTr);
  });

  el.appendChild(tbl);

  // Wire buttons via delegation
  tbl.addEventListener("click", (e) => {
    const btn = e.target.closest("button[data-action]");
    if (!btn) return;
    const idx = parseInt(btn.dataset.vidx, 10);
    const variant = result.variants.find((x) => x.idx === idx);
    if (!variant) return;

    if (btn.dataset.action === "toggle-seq") {
      const seqRow = tbl.querySelector('.gen-seq-row[data-seq-for="' + idx + '"]');
      if (!seqRow) return;
      const visible = seqRow.style.display !== "none";
      seqRow.style.display = visible ? "none" : "";
      btn.textContent = visible ? "Show seq" : "Hide seq";
    }

    if (btn.dataset.action === "send-to-screen") {
      // Fill the Screen textarea
      const ta = $("seq-input");
      if (!ta) return;
      ta.value = variant.sequence || "";
      if (typeof updateSeqMeta === "function") updateSeqMeta();

      // Switch to Screen tab
      const screenTab = document.querySelector('.tab[data-tab="screen"]');
      if (screenTab) screenTab.click();

      // Auto-click Screen after a short delay so the user sees the textarea populate
      setTimeout(() => {
        const screenBtn = $("btn-screen");
        if (screenBtn && !screenBtn.disabled) screenBtn.click();
      }, 600);
    }
  });

  // Cross-match for first variant
  if (result.variants[0] && result.variants[0].cross_ranking && result.variants[0].cross_ranking.length) {
    const cross = document.createElement("div");
    cross.className = "gen-cross";
    cross.innerHTML =
      '<div class="gen-cross-head">Cross-match: variant v0 fingerprint vs all other toxins</div>' +
      result.variants[0].cross_ranking.map(([acc, score]) => {
        const pct = Math.min(100, score * 100).toFixed(0);
        return '<div class="gen-cross-row">' +
          '<span class="nm">' + acc + '</span>' +
          '<div class="bar"><span style="width:' + pct + '%"></span></div>' +
          '<span class="sc">' + score.toFixed(4) + '</span>' +
          '</div>';
      }).join("");
    el.appendChild(cross);
  }

  // Also feed the pipeline trace tab with this generation run
  showGenerateTrace(result);
}

// Load toxins when Generate tab is clicked
document.querySelectorAll(".tab").forEach((btn) => {
  if (btn.dataset.tab === "generate") {
    btn.addEventListener("click", () => setTimeout(loadToxins, 50));
  }
});


/* ============================================================
   CONTRASTIVE PANEL — Screen tab
   ============================================================ */
function renderContrastive(j) {
  const hero = $("ctr-hero");
  if (!hero) return;

  const flagged = j.contrastive_flagged || j.flagged || false;
  hero.className = "ctr-hero " + (flagged ? "flagged" : "clear");

  $("ctr-glyph").textContent = flagged ? "\u26A0" : "\u2713";
  $("ctr-name").textContent =
    (flagged ? "Flagged" : "Clear") + " \u00B7 best match: " + (j.best_toxin_name || "—");
  $("ctr-score").textContent = (j.best_toxin_score || 0).toFixed(4);

  $("ctr-toxin").textContent =
    (j.best_toxin_name || "—") + "  " + (j.best_toxin_score || 0).toFixed(4);
  $("ctr-safe").textContent =
    (j.best_safe_name || "—") + "  " + (j.best_safe_score || 0).toFixed(4);

  const delta = j.hazard_delta || 0;
  const deltaStr = (delta >= 0 ? "+" : "") + delta.toFixed(4);
  $("ctr-delta").textContent = deltaStr;

  // Gates
  const thrDelta = j.threshold_delta || 0.02;
  const thrFloor = j.threshold_floor || 0.90;
  $("ctr-thr-delta").textContent = thrDelta.toFixed(2);
  $("ctr-thr-floor").textContent = thrFloor.toFixed(2);

  updateGateBar(
    "ctr-gate-delta",
    delta,
    thrDelta,
    (v) => (v >= 0 ? "+" : "") + v.toFixed(4)
  );
  updateGateBar(
    "ctr-gate-floor",
    j.best_toxin_score || 0,
    thrFloor,
    (v) => v.toFixed(4)
  );
}

function updateGateBar(containerId, value, threshold, format) {
  const el = $(containerId);
  if (!el) return;

  const passed = value >= threshold;
  el.classList.toggle("pass", passed);
  el.classList.toggle("fail", !passed);

  const fill = el.querySelector(".gate-fill");
  const thr  = el.querySelector(".gate-threshold");
  const val  = el.querySelector(".gate-val");

  const maxScale = Math.max(Math.abs(threshold) * 3, Math.abs(value) * 1.2, 0.01);
  const width = Math.min(100, Math.max(0, (value / maxScale) * 100));
  const thrPos = Math.min(100, Math.max(0, (threshold / maxScale) * 100));

  fill.style.width = width + "%";
  thr.style.left = thrPos + "%";
  val.textContent = format(value);
}


/* ============================================================
   CUSTOM SEQUENCE MODE — Generate tab
   ============================================================ */

// ----- Mode switch (curated vs custom) -----
document.querySelectorAll(".gen-mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    const mode = btn.dataset.genMode;
    document.querySelectorAll(".gen-mode-btn").forEach((b) =>
      b.classList.toggle("active", b === btn));
    document.querySelectorAll(".gen-subpanel").forEach((p) =>
      p.classList.toggle("active", p.dataset.genSub === mode));
  });
});

// ----- Sequence length counter -----
const customSeqEl = $("gen-custom-seq");
if (customSeqEl) {
  customSeqEl.addEventListener("input", () => {
    const s = customSeqEl.value.replace(/\s/g, "").toUpperCase();
    $("gen-custom-len").textContent = s.length + " aa";
    // Auto-set region end if not manually set
    const endEl = $("gen-custom-end");
    if (endEl && !endEl.value) {
      endEl.placeholder = "(end of seq = " + s.length + ")";
    }
  });
}

// ----- UniProt fetch -----
const fetchBtn = $("gen-uniprot-fetch");
if (fetchBtn) {
  fetchBtn.addEventListener("click", async () => {
    const acc = ($("gen-uniprot-acc").value || "").trim().toUpperCase();
    if (!acc) return;
    fetchBtn.disabled = true;
    fetchBtn.textContent = "...";
    $("gen-uniprot-meta").innerHTML = "";

    try {
      const r = await fetch("/api/uniprot_fetch/" + encodeURIComponent(acc));
      const j = await r.json();
      if (j.error) throw new Error(j.error);

      // Populate textarea
      customSeqEl.value = j.sequence;
      customSeqEl.dispatchEvent(new Event("input"));

      // Show metadata
      const siteBits = [];
      if (j.active_sites && j.active_sites.length) {
        j.active_sites.forEach((s) => siteBits.push(
          '<span class="fm-site">active: ' + s.start + (s.end !== s.start ? "–" + s.end : "") + "</span>"
        ));
      }
      if (j.binding_sites && j.binding_sites.length) {
        j.binding_sites.slice(0, 4).forEach((s) => siteBits.push(
          '<span class="fm-site">binding: ' + s.start + (s.end !== s.start ? "–" + s.end : "") + "</span>"
        ));
      }

      $("gen-uniprot-meta").innerHTML =
        '<div><span class="fm-label">accession</span><span class="fm-val">' + j.accession + "</span></div>" +
        '<div><span class="fm-label">name</span><span class="fm-val">' + (j.name || "—") + "</span></div>" +
        '<div><span class="fm-label">organism</span><span class="fm-val">' + (j.organism || "—") + "</span></div>" +
        '<div><span class="fm-label">length</span><span class="fm-val">' + j.length + " aa</span></div>" +
        (siteBits.length ? '<div style="margin-top:6px">' + siteBits.join("") + "</div>" : "");

      // If UniProt has active sites, pre-fill region
      if (j.active_sites && j.active_sites.length) {
        const first = j.active_sites[0].start;
        const last = j.active_sites[j.active_sites.length - 1].end;
        const pad = 50;
        $("gen-custom-start").value = Math.max(0, first - pad - 1);
        $("gen-custom-end").value = Math.min(j.length, last + pad);
      }
    } catch (e) {
      $("gen-uniprot-meta").innerHTML =
        '<div style="color:#f59e0b">Fetch failed: ' + e.message + "</div>";
    } finally {
      fetchBtn.disabled = false;
      fetchBtn.textContent = "Fetch";
    }
  });
}

// ----- Custom generation start -----
const customStart = $("gen-custom-start");
if (customStart) {
  customStart.addEventListener("click", async () => {
    const seq = (customSeqEl.value || "").replace(/\s/g, "").toUpperCase();
    if (seq.length < 50) {
      alert("Sequence must be at least 50 residues");
      return;
    }
    const mode = $("gen-custom-mode").value;
    const n = parseInt($("gen-custom-n").value) || 3;
    const start = parseInt($("gen-custom-start").value) || 0;
    const endRaw = $("gen-custom-end").value;
    const end = endRaw ? parseInt(endRaw) : null;
    const autoCat = $("gen-auto-catalytic").checked;

    customStart.disabled = true;
    customStart.textContent = "Generating...";

    $("gen-custom-status").style.display = "block";
    $("gen-custom-stage").textContent = "queued";
    $("gen-custom-pct").textContent = "0%";
    $("gen-custom-bar").style.width = "0%";
    $("gen-custom-note").textContent = mode === "aggressive"
      ? "Aggressive mode redesigns the entire region."
      : "Preserved mode protects auto-detected catalytic-like residues.";
    $("gen-custom-results").innerHTML = "";

    try {
      const r = await fetch("/api/generate_custom", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence: seq,
          mode: mode,
          n_variants: n,
          region_start: start,
          region_end: end,
          auto_catalytic: autoCat,
        }),
      });
      const j = await r.json();
      if (j.error) throw new Error(j.error);
      pollCustomJob(j.job_id);
    } catch (e) {
      $("gen-custom-stage").textContent = "error: " + e.message;
      customStart.disabled = false;
      customStart.textContent = "Generate live";
    }
  });
}

let _customPollTimer = null;
function pollCustomJob(jobId) {
  if (_customPollTimer) clearInterval(_customPollTimer);
  _customPollTimer = setInterval(async () => {
    try {
      const r = await fetch("/api/generate/" + jobId);
      const job = await r.json();

      if (job.status === "ready") {
        clearInterval(_customPollTimer);
        $("gen-custom-stage").textContent = "done in " + (job.elapsed_s || "?") + " s";
        $("gen-custom-pct").textContent = "100%";
        $("gen-custom-bar").style.width = "100%";
        $("gen-custom-start").disabled = false;
        $("gen-custom-start").textContent = "Generate live";
        renderGenerateResult(job.result, "gen-custom-results");
      } else if (job.status === "error") {
        clearInterval(_customPollTimer);
        $("gen-custom-stage").textContent = "error: " + (job.error || "unknown");
        $("gen-custom-start").disabled = false;
        $("gen-custom-start").textContent = "Generate live";
      } else {
        $("gen-custom-stage").textContent = job.stage || job.status;
        const pct = Math.round((job.progress || 0) * 100);
        $("gen-custom-pct").textContent = pct + "%";
        $("gen-custom-bar").style.width = pct + "%";
      }
    } catch (e) {
      console.error("poll error:", e);
    }
  }, 2000);
}

/* ============================================================
   REFACTOR: make renderGenerateResult target-agnostic
   (Accept an optional container ID; default to gen-results)
   ============================================================ */
if (typeof window._origRenderGenerateResult === "undefined") {
  window._origRenderGenerateResult = window.renderGenerateResult;
}

window.renderGenerateResult = function(result, targetId) {
  const el = $(targetId || "gen-results");
  if (!el) return;
  el.innerHTML = "";

  // Header
  const head = document.createElement("div");
  head.className = "gen-cross";
  head.innerHTML =
    '<div class="gen-cross-head">' +
      result.name + ' &middot; ' + result.mode + ' mode &middot; ' +
      result.region_length + ' aa region &middot; ' +
      result.n_catalytic + ' catalytic residues' +
      (result.source === "custom" ? ' &middot; custom sequence' : '') +
      (result.auto_catalytic ? ' &middot; auto-detected' : '') +
    '</div>';
  el.appendChild(head);

  // Results table (same layout as curated)
  const tbl = document.createElement("table");
  tbl.className = "gen-result-table";
  tbl.innerHTML =
    '<thead><tr>' +
      '<th>Variant</th>' +
      '<th>Identity</th>' +
      '<th>BLAST</th>' +
      '<th>Contrastive</th>' +
      '<th>Fingerprint</th>' +
      '<th>Verdict</th>' +
      '<th>Signal</th>' +
      '<th>Tier / Premium</th>' +
      '<th>Actions</th>' +
    '</tr></thead><tbody></tbody>';

  const tbody = tbl.querySelector("tbody");
  result.variants.forEach((v) => {
    const tr = document.createElement("tr");
    tr.className = "gen-variant-row";

    const blCls = v.blast_evaded ? "evaded" : "blast-caught";
    const fpCls = v.fingerprint_caught ? "caught" : "clear";
    const ctrFlagged = v.contrastive_flagged || false;
    const ctrCls = ctrFlagged ? "caught" : "clear";
    const ctrLabel = ctrFlagged ? "FLAG" : "clear";

    const verdict = v.verdict || "—";
    const verdictCls = verdict === "FLAGGED" ? "evaded" : "caught";
    const tier = (v.pricing && v.pricing.tier) ? v.pricing.tier : "—";
    const premium = (v.pricing && v.pricing.premium_usd !== undefined)
      ? "$" + v.pricing.premium_usd.toFixed(2) : "—";

    const ctrDelta = (v.contrastive_delta !== undefined && v.contrastive_delta !== null)
      ? (v.contrastive_delta >= 0 ? "+" : "") + v.contrastive_delta.toFixed(4)
      : "—";

    const reason = v.verdict_reason || "—";
    let sigClass = "clear", sigLabel = "—";
    if (reason.indexOf("Both") === 0) { sigClass = "caught"; sigLabel = "BOTH"; }
    else if (reason.indexOf("Fingerprint") === 0) { sigClass = "caught"; sigLabel = "FP"; }
    else if (reason.indexOf("Contrastive") === 0) { sigClass = "caught"; sigLabel = "CTR"; }
    else if (reason.indexOf("Neither") === 0) { sigClass = "clear"; sigLabel = "NONE"; }

    const blastCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + v.blast.toFixed(3) + '</span>' +
        '<span class="gen-badge ' + blCls + '">' + (v.blast_evaded ? 'EVADED' : 'caught') + '</span>' +
      '</div>';
    const ctrCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + ctrDelta + '</span>' +
        '<span class="gen-badge ' + ctrCls + '">' + ctrLabel + '</span>' +
      '</div>';
    const fpCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + (v.fingerprint !== null ? v.fingerprint.toFixed(3) : '—') + '</span>' +
        '<span class="gen-badge ' + fpCls + '">' + (v.fingerprint_caught ? 'CAUGHT' : 'clear') + '</span>' +
      '</div>';
    const tierCell =
      '<div class="gen-stack">' +
        '<span class="gen-num">' + tier + '</span>' +
        '<span class="gen-num-muted">' + premium + '</span>' +
      '</div>';

    tr.innerHTML =
      '<td>v' + v.idx + '</td>' +
      '<td>' + (v.identity * 100).toFixed(2) + '%</td>' +
      '<td>' + blastCell + '</td>' +
      '<td>' + ctrCell + '</td>' +
      '<td>' + fpCell + '</td>' +
      '<td><span class="gen-badge ' + verdictCls + '">' + verdict + '</span></td>' +
      '<td><span class="gen-badge ' + sigClass + '">' + sigLabel + '</span></td>' +
      '<td>' + tierCell + '</td>' +
      '<td>' +
        '<button class="gen-mini" data-vidx="' + v.idx + '" data-action="toggle-seq">seq</button> ' +
        '<button class="gen-mini primary" data-vidx="' + v.idx + '" data-action="send-to-screen">→ Screen</button>' +
      '</td>';
    tbody.appendChild(tr);

    // Hidden sequence row
    const seqTr = document.createElement("tr");
    seqTr.className = "gen-seq-row";
    seqTr.dataset.seqFor = v.idx;
    seqTr.style.display = "none";
    seqTr.innerHTML = '<td colspan="9"><pre class="gen-seq-block">' +
      (v.sequence || "") + '</pre></td>';
    tbody.appendChild(seqTr);
  });

  el.appendChild(tbl);

  // Wire buttons
  tbl.addEventListener("click", (e) => {
    const btn = e.target.closest("button[data-action]");
    if (!btn) return;
    const idx = parseInt(btn.dataset.vidx, 10);
    const variant = result.variants.find((x) => x.idx === idx);
    if (!variant) return;

    if (btn.dataset.action === "toggle-seq") {
      const seqRow = tbl.querySelector('.gen-seq-row[data-seq-for="' + idx + '"]');
      if (!seqRow) return;
      const visible = seqRow.style.display !== "none";
      seqRow.style.display = visible ? "none" : "";
      btn.textContent = visible ? "seq" : "hide";
    }

    if (btn.dataset.action === "send-to-screen") {
      const ta = $("seq-input");
      if (!ta) return;
      ta.value = variant.sequence || "";
      if (typeof updateSeqMeta === "function") updateSeqMeta();
      const screenTab = document.querySelector('.tab[data-tab="screen"]');
      if (screenTab) screenTab.click();
      setTimeout(() => {
        const screenBtn = $("btn-screen");
        if (screenBtn && !screenBtn.disabled) screenBtn.click();
      }, 600);
    }
  });

  // Cross-match for first variant
  if (result.variants[0] && result.variants[0].cross_ranking && result.variants[0].cross_ranking.length) {
    const cross = document.createElement("div");
    cross.className = "gen-cross";
    cross.innerHTML =
      '<div class="gen-cross-head">Cross-match: variant v0 fingerprint vs all other toxins</div>' +
      result.variants[0].cross_ranking.map(([acc, score]) => {
        const pct = Math.min(100, score * 100).toFixed(0);
        return '<div class="gen-cross-row">' +
          '<span class="nm">' + acc + '</span>' +
          '<div class="bar"><span style="width:' + pct + '%"></span></div>' +
          '<span class="sc">' + score.toFixed(4) + '</span>' +
          '</div>';
      }).join("");
    el.appendChild(cross);
  }

  // Pipeline trace
  if (result.trace) {
    showGenerateTrace(result);
  }
};


/* ============================================================
   SCREEN TAB — reference dropdowns
   ============================================================ */
let _referencesLoaded = false;

async function loadReferences() {
  if (_referencesLoaded) return;
  const tsel = $("ref-toxin");
  const ssel = $("ref-safe");
  const btn  = $("ref-load");
  if (!tsel || !ssel || !btn) return;

  try {
    const r = await fetch("/api/references");
    if (!r.ok) throw new Error("status " + r.status);
    const j = await r.json();

    // Clear and populate toxin dropdown
    tsel.innerHTML = '<option value="">select a toxin…</option>';
    j.toxins.forEach((t) => {
      const opt = document.createElement("option");
      opt.value = t.accession;
      opt.textContent = t.name + "  (" + t.accession + ")";
      tsel.appendChild(opt);
    });

    // Clear and populate safe dropdown
    ssel.innerHTML = '<option value="">select a safe protein…</option>';
    j.safes.forEach((t) => {
      const opt = document.createElement("option");
      opt.value = t.accession;
      opt.textContent = t.name + "  (" + t.accession + ")";
      ssel.appendChild(opt);
    });

    _referencesLoaded = true;
  } catch (e) {
    console.error("references fetch failed:", e);
  }
}

// Enable load button when either dropdown has a value
function updateRefLoadBtn() {
  const tsel = $("ref-toxin");
  const ssel = $("ref-safe");
  const btn  = $("ref-load");
  if (!btn) return;
  const hasToxin = tsel && tsel.value;
  const hasSafe  = ssel && ssel.value;
  // Only one at a time — prefer whichever was changed most recently
  btn.disabled = !(hasToxin || hasSafe);
}

// When user changes one dropdown, clear the other
const tsel = $("ref-toxin");
const ssel = $("ref-safe");
const loadBtn = $("ref-load");

if (tsel) {
  tsel.addEventListener("change", () => {
    if (tsel.value && ssel) ssel.value = "";
    updateRefLoadBtn();
  });
}
if (ssel) {
  ssel.addEventListener("change", () => {
    if (ssel.value && tsel) tsel.value = "";
    updateRefLoadBtn();
  });
}

if (loadBtn) {
  loadBtn.addEventListener("click", async () => {
    let acc = (tsel && tsel.value) || (ssel && ssel.value);
    if (!acc) return;

    const originalText = loadBtn.textContent;
    loadBtn.disabled = true;
    loadBtn.textContent = "…";

    try {
      const r = await fetch("/api/reference_sequence/" + encodeURIComponent(acc));
      const j = await r.json();
      if (j.error) throw new Error(j.error);

      const ta = $("seq-input");
      if (!ta) return;
      ta.value = j.sequence || "";
      if (typeof updateSeqMeta === "function") updateSeqMeta();

      // Auto-run Screen after a short delay
      setTimeout(() => {
        const screenBtn = $("btn-screen");
        if (screenBtn && !screenBtn.disabled) screenBtn.click();
      }, 400);
    } catch (e) {
      alert("Could not load reference: " + e.message);
    } finally {
      loadBtn.textContent = originalText;
      updateRefLoadBtn();
    }
  });
}

// Preload the reference list on page load (async — no need to block)
setTimeout(loadReferences, 800);
