/* Minimal smoke + numerical-sanity test for the KRAS G12C Living Meta dashboard.
 *
 * The statistical engine ships inside the single-file HTML app (it runs in the
 * browser inside IIFEs and is not importable), so this test does two things:
 *
 *   1. Structural integrity of the shipped assets (no BOM, redirect target
 *      present, no literal </script> inside JS template literals, no leaked
 *      asthma/Type-2-inflammation template strings, balanced enough markup).
 *   2. Numerical sanity: an independent re-implementation of the exact pooling
 *      math the app uses (continuity-corrected log-OR, inverse-variance fixed
 *      pool, HR-from-published-CI, Paule-Mandel tau^2) checked against
 *      hand-computed reference values. If the app's math is changed in an
 *      incompatible way, the reference values here document the contract.
 *
 * Run: node test/smoke.test.js   (exit 0 = all pass, exit 1 = failure)
 */
'use strict';
const fs = require('fs');
const path = require('path');

let passed = 0, failed = 0;
function ok(name, cond, extra) {
  if (cond) { passed++; console.log('  ok  - ' + name); }
  else { failed++; console.log('  FAIL- ' + name + (extra ? '  (' + extra + ')' : '')); }
}
function approx(a, b, tol) { return Math.abs(a - b) <= (tol == null ? 1e-6 : tol); }

const root = path.join(__dirname, '..');
const mainHtml = fs.readFileSync(path.join(root, 'KRAS_G12C_NSCLC_REVIEW.html'), 'utf8');
const indexHtml = fs.readFileSync(path.join(root, 'index.html'), 'utf8');

console.log('# structural integrity');
// No UTF-8 BOM in shipped assets.
ok('main html has no BOM', mainHtml.charCodeAt(0) !== 0xFEFF);
ok('index html has no BOM', indexHtml.charCodeAt(0) !== 0xFEFF);
// Redirect target exists.
ok('index redirects to the dashboard file',
   indexHtml.indexOf('KRAS_G12C_NSCLC_REVIEW.html') !== -1);
// No leaked asthma / Type-2-inflammation template strings (wrong domain).
ok('no leaked "Type 2 Inflammation" string', mainHtml.indexOf('Type 2 Inflammation') === -1);
ok('no leaked "IL-4Ralpha" string', mainHtml.indexOf('IL-4Ralpha') === -1);
ok('no leaked "RapidMeta Respiratory" title', mainHtml.indexOf('RapidMeta Respiratory') === -1);
ok('no leaked eosinophil subgroup value', mainHtml.indexOf('Blood eosinophils, smoking status') === -1);
ok('no leaked "type-2 inflammation" descriptor', /type[- ]?2 inflammation/i.test(mainHtml) === false);
ok('no leaked "exacerbation" (asthma) endpoint', /exacerbation/i.test(mainHtml) === false);
ok('no doubled "Adults with Adults" template artifact', mainHtml.indexOf('Adults with Adults') === -1);
ok('no leaked "RapidMeta Cardiology" title', mainHtml.indexOf('RapidMeta Cardiology') === -1);
// ClinicalTrials.gov live query is clean (no leaked PICO fragment in intr).
ok('ctgov intervention query is clean',
   mainHtml.indexOf('query.intr=kras g12c AND Adults') === -1);
ok('ctgov intervention query targets KRAS G12C',
   mainHtml.indexOf('query.intr=KRAS G12C') !== -1);
// <script> open/close tags are balanced (corruption / unterminated-block guard).
{
  const openTags = (mainHtml.match(/<script[\s>]/g) || []).length;
  const closeTags = (mainHtml.match(/<\/script>/g) || []).length;
  ok('<script> open/close tags are balanced', openTags === closeTags,
     openTags + ' open vs ' + closeTags + ' close');
}
// No unfilled placeholder tokens in shipped markup.
// (token strings are assembled so this negative-test fixture does not itself
//  trip placeholder linters)
{
  const tok1 = 'REPLACE' + '_ME';
  const tok2 = '__PLACE' + 'HOLDER__';
  ok('no unfilled placeholder tokens',
     mainHtml.indexOf(tok1) === -1 && mainHtml.indexOf(tok2) === -1);
}

console.log('# numerical sanity (independent re-implementation of the app math)');

// --- continuity-corrected log-OR (matches binaryEffect, measure='OR') ---
function logOR(tE, tN, cE, cN) {
  const tNonE = tN - tE, cNonE = cN - cE;
  const add = (tE === 0 || tNonE === 0 || cE === 0 || cNonE === 0) ? 0.5 : 0;
  const a = tE + add, b = tNonE + add, c = cE + add, e = cNonE + add;
  return { y: Math.log(a * e / (b * c)), v: 1 / a + 1 / b + 1 / c + 1 / e };
}
// 2x2: tx 10/100, ctrl 20/100. No zero cell -> no correction.
let r = logOR(10, 100, 20, 100);
// OR = (10*80)/(90*20) = 800/1800 = 0.4444...; logOR = ln(0.44444)
ok('logOR point value (no zero cell)', approx(r.y, Math.log((10 * 80) / (90 * 20))));
ok('logOR variance (no zero cell)', approx(r.v, 1 / 10 + 1 / 90 + 1 / 20 + 1 / 80));
// zero-cell triggers 0.5 correction.
let rz = logOR(0, 50, 5, 50);
ok('zero cell applies 0.5 continuity correction',
   approx(rz.y, Math.log((0.5 * 45.5) / (50.5 * 5.5))));

// --- HR from published CI (matches hrFromPublished) ---
function hrFromCI(hr, lo, hi) {
  return { y: Math.log(hr), v: Math.pow((Math.log(hi) - Math.log(lo)) / (2 * 1.959964), 2) };
}
// CodeBreaK 200 PFS HR 0.66 (e.g. 95% CI 0.51-0.86) -> log-scale SE check.
let h = hrFromCI(0.66, 0.51, 0.86);
ok('logHR point = ln(HR)', approx(h.y, Math.log(0.66)));
ok('SE-from-CI is positive and finite', h.v > 0 && isFinite(h.v));

// --- inverse-variance fixed pool (matches fePool / poolWith with tau2=0) ---
function fePool(ys, vs) {
  const w = vs.map(v => 1 / v);
  const sW = w.reduce((s, wi) => s + wi, 0);
  const mu = ys.reduce((s, y, i) => s + y * w[i], 0) / sW;
  return { mu, se: Math.sqrt(1 / sW) };
}
// Two identical effects -> pooled equals the common effect; SE shrinks by sqrt(2).
let p = fePool([Math.log(0.6), Math.log(0.6)], [0.04, 0.04]);
ok('FE pool of identical effects equals common effect', approx(p.mu, Math.log(0.6)));
ok('FE pool SE = sqrt(v/2) for two equal-variance studies',
   approx(p.se, Math.sqrt(0.04 / 2)));

// --- Paule-Mandel tau^2 (matches pauleMandelTau2) ---
function pmTau2(y, v) {
  const k = y.length;
  if (k < 2) return 0;
  let tau2 = 0;
  const target = k - 1;
  for (let it = 0; it < 200; it++) {
    const w = v.map(vi => 1 / (vi + tau2));
    const sW = w.reduce((a, b) => a + b, 0);
    let yBar = 0; for (let i = 0; i < k; i++) yBar += y[i] * w[i]; yBar /= sW;
    let Q = 0; for (let j = 0; j < k; j++) { const d = y[j] - yBar; Q += w[j] * d * d; }
    const diff = Q - target;
    if (Math.abs(diff) < 1e-7) break;
    let slope = 0; for (let m = 0; m < k; m++) { const dm = y[m] - yBar; slope += -w[m] * w[m] * dm * dm; }
    if (!isFinite(slope) || Math.abs(slope) < 1e-14) break;
    const next = Math.max(0, tau2 - diff / slope);
    if (Math.abs(next - tau2) < 1e-10) { tau2 = next; break; }
    tau2 = next;
  }
  return tau2;
}
// Homogeneous studies (Q <= k-1) -> tau^2 clamped to 0.
ok('PM tau^2 = 0 for homogeneous studies',
   approx(pmTau2([0.0, 0.0, 0.0], [0.05, 0.05, 0.05]), 0, 1e-9));
// Heterogeneous studies -> tau^2 strictly positive.
let t2 = pmTau2([Math.log(0.4), Math.log(1.6), Math.log(0.5)], [0.02, 0.02, 0.02]);
ok('PM tau^2 > 0 for heterogeneous studies', t2 > 0, 'tau2=' + t2);
// k=1 guard.
ok('PM tau^2 = 0 when k<2 (guard)', pmTau2([0.3], [0.05]) === 0);

console.log('\n' + passed + ' passed, ' + failed + ' failed');
if (failed > 0) process.exit(1);
