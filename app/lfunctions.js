/* lfunctions.js — evaluation of the example L-functions in the critical
 * strip, in double precision, fast enough to colour a grid.
 *
 * Normalisation is arithmetic throughout (the one LMFDB calls "arithmetic"):
 *   zeta, Dirichlet L and Dedekind zeta have their critical line at Re(s)=1/2;
 *   the elliptic-curve L-functions have theirs at Re(s)=1, with the functional
 *   equation s <-> 2-s.
 *
 * Methods:
 *   zeta, L(s,chi), zeta_K   Euler-Maclaurin on the Hurwitz zeta function.
 *                            Valid in the whole plane (bar the pole at s=1)
 *                            because the formula IS the analytic continuation,
 *                            not a rearranged series.
 *   L(E,s)                   the smoothed sum over the incomplete gamma
 *                            function, Lambda(s) = sum a_n [x^-s G(s,x)
 *                            + eps x^(s-2) G(2-s,x)], x = 2 pi n / sqrt(N).
 *
 * Every routine here is checked against Magma 2.29-10 in test-lfunctions.html.
 */
(function (global) {
  'use strict';

  var TWO_PI = 2 * Math.PI;

  /* ---------------------------------------------------------- complex ---- */

  function cmul(ar, ai, br, bi) { return [ar * br - ai * bi, ar * bi + ai * br]; }

  function cdiv(ar, ai, br, bi) {
    var d = br * br + bi * bi;
    return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d];
  }

  function clog(ar, ai) { return [Math.log(Math.hypot(ar, ai)), Math.atan2(ai, ar)]; }

  function cexp(ar, ai) { var e = Math.exp(ar); return [e * Math.cos(ai), e * Math.sin(ai)]; }

  function csin(ar, ai) {
    return [Math.sin(ar) * Math.cosh(ai), Math.cos(ar) * Math.sinh(ai)];
  }

  /* x^z for real x > 0 */
  function rpow(x, zr, zi) {
    var L = Math.log(x), e = Math.exp(zr * L), th = zi * L;
    return [e * Math.cos(th), e * Math.sin(th)];
  }

  /* ------------------------------------------------------------ gamma ---- */

  var LANCZOS = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  ];
  var SQRT_2PI = Math.sqrt(TWO_PI);

  function cgamma(zr, zi) {
    if (zr < 0.5) {
      var s = csin(Math.PI * zr, Math.PI * zi);
      var g = cgamma(1 - zr, -zi);
      var den = cmul(s[0], s[1], g[0], g[1]);
      return cdiv(Math.PI, 0, den[0], den[1]);
    }
    var wr = zr - 1, wi = zi;
    var xr = LANCZOS[0], xi = 0;
    for (var i = 1; i < 9; i++) {
      var t = cdiv(LANCZOS[i], 0, wr + i, wi);
      xr += t[0]; xi += t[1];
    }
    var tr = wr + 7.5, ti = wi;
    var lt = clog(tr, ti);
    var ex = cmul(wr + 0.5, wi, lt[0], lt[1]);
    var p = cexp(ex[0] - tr, ex[1] - ti);
    var r = cmul(p[0], p[1], xr, xi);
    return [SQRT_2PI * r[0], SQRT_2PI * r[1]];
  }

  /* -------------------------------------------- upper incomplete gamma --- */
  /* G(a, x) for complex a and real x > 0.
   * Legendre's continued fraction when x dominates |a| (Lentz's algorithm),
   * otherwise the lower-series plus G(a). The switch also keeps us off the
   * poles of G(a): those only matter in the series branch, and that branch is
   * never taken when |a| is small, since then x > |a| + 1 holds. */
  function incGammaUpper(ar, ai, x) {
    var amod = Math.hypot(ar, ai);
    var TINY = 1e-300, i, del;

    if (x > amod + 1) {
      var br = x + 1 - ar, bi = -ai;
      var cr = 1 / TINY, ci = 0;
      var d = cdiv(1, 0, br, bi);
      var hr = d[0], hi = d[1];
      for (i = 1; i <= 400; i++) {
        var anr = -i * (i - ar), ani = i * ai;
        br += 2;
        var ad = cmul(anr, ani, d[0], d[1]);
        var dr = br + ad[0], di = bi + ad[1];
        if (Math.abs(dr) + Math.abs(di) < TINY) { dr = TINY; di = 0; }
        d = cdiv(1, 0, dr, di);
        var ac = cdiv(anr, ani, cr, ci);
        cr = br + ac[0]; ci = bi + ac[1];
        if (Math.abs(cr) + Math.abs(ci) < TINY) { cr = TINY; ci = 0; }
        del = cmul(d[0], d[1], cr, ci);
        var nh = cmul(hr, hi, del[0], del[1]);
        hr = nh[0]; hi = nh[1];
        if (Math.abs(del[0] - 1) + Math.abs(del[1]) < 1e-16) break;
      }
      var pf = rpow(x, ar, ai), ex = Math.exp(-x);
      return [ex * (pf[0] * hr - pf[1] * hi), ex * (pf[0] * hi + pf[1] * hr)];
    }

    var apr = ar, api = ai;
    var sum = cdiv(1, 0, ar, ai);
    var dr2 = sum[0], di2 = sum[1];
    for (i = 1; i <= 800; i++) {
      apr += 1;
      var t = cdiv(dr2 * x, di2 * x, apr, api);
      dr2 = t[0]; di2 = t[1];
      sum[0] += dr2; sum[1] += di2;
      if (Math.hypot(dr2, di2) < Math.hypot(sum[0], sum[1]) * 1e-17) break;
    }
    var pf2 = rpow(x, ar, ai), e2 = Math.exp(-x);
    var lowR = e2 * (pf2[0] * sum[0] - pf2[1] * sum[1]);
    var lowI = e2 * (pf2[0] * sum[1] + pf2[1] * sum[0]);
    var G = cgamma(ar, ai);
    return [G[0] - lowR, G[1] - lowI];
  }

  /* ---------------------------------------------------- Hurwitz zeta ----- */
  /* B_{2j} / (2j)! for j = 1..10 */
  var EM_C = [
    1 / 12, -1 / 720, 1 / 30240, -1 / 1209600, 5 / 239500800,
    -691 / 1307674368000, 7 / 523069747200, -3617 / 10670622842880000,
    43867 / 5109094217170944000, -174611 / 802857662698291200000
  ];

  function emTerms(si) { return Math.max(16, Math.ceil(Math.abs(si)) + 16); }

  /* sum_k c_k Q_k^{1-s} / (s-1).
   * When sum c_k = 0 -- which is exactly the case for a non-principal
   * character -- this is finite at s = 1 even though every term blows up, so
   * near the pole the cancellation is done symbolically through the series
   * rather than numerically. */
  function poleCombo(sr, si, Qs, cs) {
    var wr = sr - 1, wi = si, k;
    if (Math.hypot(wr, wi) > 0.05) {
      var re = 0, im = 0;
      for (k = 0; k < Qs.length; k++) {
        var t = rpow(Qs[k], 1 - sr, -si);
        re += cs[k] * t[0]; im += cs[k] * t[1];
      }
      return cdiv(re, im, wr, wi);
    }
    var L = [], pw = [];
    for (k = 0; k < Qs.length; k++) { L.push(Math.log(Qs[k])); pw.push(1); }
    var wpr = 1, wpi = 0, fact = 1, R = 0, I = 0, sgn = -1;
    for (var j = 1; j <= 24; j++) {
      fact *= j;
      var S = 0;
      for (k = 0; k < Qs.length; k++) { pw[k] *= L[k]; S += cs[k] * pw[k]; }
      var c = sgn * S / fact;
      R += c * wpr; I += c * wpi;
      var nr = wpr * wr - wpi * wi;
      wpi = wpr * wi + wpi * wr; wpr = nr;
      sgn = -sgn;
    }
    return [R, I];
  }

  /* Everything in Euler-Maclaurin for zeta(s,a) EXCEPT the pole term. */
  function hurwitzCore(sr, si, a, N) {
    var re = 0, im = 0, k, t;
    for (k = 0; k < N; k++) {
      t = rpow(k + a, -sr, -si);
      re += t[0]; im += t[1];
    }
    var q = N + a;
    var w = rpow(q, -sr, -si);
    re += 0.5 * w[0];
    im += 0.5 * w[1];
    var pr = sr, pi = si, qp = 1 / q, qq = q * q;
    for (var j = 1; j <= EM_C.length; j++) {
      var term = cmul(pr, pi, w[0], w[1]);
      re += EM_C[j - 1] * term[0] * qp;
      im += EM_C[j - 1] * term[1] * qp;
      var m = cmul(pr, pi, sr + 2 * j - 1, si);
      m = cmul(m[0], m[1], sr + 2 * j, si);
      pr = m[0]; pi = m[1];
      qp /= qq;
    }
    return [re, im];
  }

  function hurwitzZeta(sr, si, a, N) {
    if (N === undefined) N = emTerms(si);
    var c = hurwitzCore(sr, si, a, N);
    var w = rpow(N + a, 1 - sr, -si);
    var p = cdiv(w[0], w[1], sr - 1, si);
    return [c[0] + p[0], c[1] + p[1]];
  }

  function zeta(sr, si) { return hurwitzZeta(sr, si, 1); }

  /* ------------------------------------------------ Dirichlet L-series --- */
  /* L(s,chi) = q^{-s} sum_{a=1..q} chi(a) zeta(s, a/q).
   *
   * Written in terms of Q_a = qN + a, every q^{-s} folds into the Hurwitz
   * terms and the whole thing becomes
   *
   *   sum_{n<=qN} chi(n) n^{-s}                      (an ordinary partial sum)
   * + (1/q) sum_a chi(a) Q_a^{1-s}/(s-1)             (poles, which cancel)
   * + sum_a chi(a) [ Q_a^{-s}/2 + Bernoulli terms ]
   *
   * so only powers n^{-s} are ever needed, which is what makes the grid
   * version below fast. */
  function dirichletL(sr, si, ch, N) {
    if (N === undefined) N = emTerms(si);
    var q = ch.q, v = ch.v, M = q * N, n, t, a;
    var re = 0, im = 0;
    for (n = 1; n <= M; n++) {
      var c = v[(n - 1) % q];
      if (!c) continue;
      t = rpow(n, -sr, -si);
      re += c * t[0]; im += c * t[1];
    }
    var Qs = [], cs = [];
    for (a = 1; a <= q; a++) {
      if (!v[a - 1]) continue;
      Qs.push(M + a); cs.push(v[a - 1] / q);
    }
    var P = poleCombo(sr, si, Qs, cs);
    re += P[0]; im += P[1];

    for (a = 1; a <= q; a++) {
      var ca = v[a - 1];
      if (!ca) continue;
      var Q = M + a;
      var w = rpow(Q, -sr, -si);
      re += ca * 0.5 * w[0];
      im += ca * 0.5 * w[1];
      var pr = sr, pi = si, f = q / Q, fp = f;
      for (var j = 1; j <= EM_C.length; j++) {
        var term = cmul(pr, pi, w[0], w[1]);
        re += ca * EM_C[j - 1] * term[0] * fp;
        im += ca * EM_C[j - 1] * term[1] * fp;
        var m = cmul(pr, pi, sr + 2 * j - 1, si);
        m = cmul(m[0], m[1], sr + 2 * j, si);
        pr = m[0]; pi = m[1];
        fp *= f * f;
      }
    }
    return [re, im];
  }

  /* For a quadratic field of discriminant D, zeta_K = zeta * L(s, chi_D). */
  function dedekind(sr, si, ch) {
    var z = zeta(sr, si), l = dirichletL(sr, si, ch);
    return cmul(z[0], z[1], l[0], l[1]);
  }

  /* ------------------------------------------- elliptic curve L-series --- */

  function ellipticL(sr, si, curve) {
    var Nc = curve.N, eps = curve.sign, a = curve.a;
    var c = TWO_PI / Math.sqrt(Nc);
    var T = Math.min(a.length - 1, curve.terms);
    var Lr = 0, Li = 0;
    for (var n = 1; n <= T; n++) {
      var an = a[n];
      if (an === 0) continue;
      var x = c * n;
      var g1 = incGammaUpper(sr, si, x);
      var p1 = rpow(x, -sr, -si);
      var g2 = incGammaUpper(2 - sr, -si, x);
      var p2 = rpow(x, sr - 2, si);
      var t1 = cmul(p1[0], p1[1], g1[0], g1[1]);
      var t2 = cmul(p2[0], p2[1], g2[0], g2[1]);
      Lr += an * (t1[0] + eps * t2[0]);
      Li += an * (t1[1] + eps * t2[1]);
    }
    var xp = rpow(c, sr, si);
    var num = cmul(Lr, Li, xp[0], xp[1]);
    var G = cgamma(sr, si);
    return cdiv(num[0], num[1], G[0], G[1]);
  }

  /* --------------------------------------------------------- the data --- */

  /* Quadratic characters chi_D, values on 1..q. From Magma's
   * KroneckerCharacter(D); all three are odd, all have root number +1. */
  var CHARS = {
    chi3:  { q: 3,  D: -3,  h: 1, w: 6,
             v: [1, -1, 0] },
    chi20: { q: 20, D: -20, h: 2, w: 2,
             v: [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, -1, 0, -1, 0, 0, 0, -1, 0, -1, 0] },
    chi23: { q: 23, D: -23, h: 3, w: 2,
             v: [1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1,
                 -1, -1, -1, -1, 0] }
  };

  var CURVES = {
    '11a1': {
      N: 11, sign: 1, terms: 30,
      a: [0, 1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0,
          2, 2, -2, -1, 0, -4, -8, 5, -4, 0, 2, 7, 8, -1, 4, -2, -4, 3, 0, -4,
          0, -8, -4, -6, 2, -2, 2, 8, 4, -3, 8, 2, 8, -6, -10, 1, 0, 0, 0, 5,
          -2, 12, -14, 4, -8, 4, 2, -7, -4, 1, 4]
    },
    '37a1': {
      N: 37, sign: -1, terms: 55,
      a: [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0,
          -4, 3, 10, 2, 0, -1, 4, -9, -2, 6, -12, -4, 8, 15, 0, 2, 12, -1, 0, 6,
          0, -9, -6, 2, -10, -12, -4, -9, 12, -6, 2, 0, -4, 1, 18, 10, 0, 0,
          -12, 8, 12, -8, 8, -6, -8, 4, -30, 8, 0, -6, -4]
    }
  };

  /* Imaginary parts of zeros on the critical line, positive ones only.
   *   zeta          Odlyzko's zeros2 table (1000-digit values, truncated)
   *   chi3          Magma 2.29-10, matching LMFDB 1-3-3.2-r1-0-0 to 3 dp
   *   chi20         Magma, matching LMFDB 1-20-20.19-r1-0-0 to all 27 digits
   *   chi23         Magma, matching LMFDB 1-23-23.22-r1-0-0 to all 27 digits
   *   11a1, 37a1    LMFDB 2-11-1.1-c1-0-0 and 2-37-1.1-c1-0-1
   * A Dedekind zeta of a quadratic field factors as zeta * L(s,chi_D), so its
   * zeros are the union of the two lists; that is assembled below. */
  var ZEROS = {
    zeta: [14.134725141734694, 21.022039638771555, 25.010857580145688,
           30.424876125859513, 32.935061587739190, 37.586178158825671,
           40.918719012147495, 43.327073280914999, 48.005150881167160,
           49.773832477672302, 52.970321477714460, 56.446247697063395,
           59.347044002602353, 60.831778524609810, 65.112544048081607,
           67.079810529494174, 69.546401711173979, 72.067157674481908,
           75.704690699083933, 77.144840068874805],
    chi3: [8.039737155681467, 11.249206207772935, 15.704619176721625,
           18.261997495693127, 20.455770807742492, 24.059414856493450,
           26.577868735774585, 28.218164506233386, 30.745040261382495,
           33.897388927259419, 35.608412653938634, 37.551796556364627,
           39.485207260929350, 42.616379226157567, 44.120572912072202,
           46.274118023513140, 47.514104510117322],
    chi20: [2.358934994086656, 4.675507749842080, 7.429109774584178,
            8.804527424544900, 10.663300734953622, 12.802400810866629,
            14.336169007500515, 15.493534961243027, 17.600734791270634,
            19.120252802798512, 20.542108704454735, 21.622206406481371,
            23.786167100158956, 24.906610399535413, 26.378159592740234,
            27.429410778006436, 29.315636462947542, 30.813432186697565,
            31.629850769137192, 33.162006756974157],
    chi23: [2.871339848930368, 4.215189804229719, 6.731189150719542,
            8.334849030124398, 10.633871230218586, 12.581696786969625,
            13.604131724738555, 15.385560733598581, 15.938101367347983,
            18.985931202695874, 19.987766847887368, 21.203602256396255,
            22.874076023463120, 23.946447295267818, 25.462561010482498,
            26.486391237836839, 28.502865498702425, 30.101719258358329,
            31.254693951229937, 31.918338693736597],
    '11a1': [6.362613894713089, 8.603539619290756, 10.035509097181079,
             11.451258610345211, 13.568639057129995, 15.914072603300384,
             17.033610320380624, 17.941433573459341, 19.185724971852241,
             20.379260464350109, 22.172490291474671, 23.301415502229715,
             25.209868424251594, 25.876403079348791, 27.067635233438503],
    '37a1': [0, 5.003170014006659, 6.870391216954432, 8.014330807872879,
             9.933098353605352, 10.775138162540800, 11.757324722849776,
             12.958386413882846, 15.603857873204318, 16.192017416874481,
             17.141693648014874, 18.063654202910710, 18.787195624663916,
             19.814822245363376]
  };
  function bySize(p, q) { return p - q; }
  ZEROS.zetaK3  = ZEROS.zeta.concat(ZEROS.chi3).sort(bySize);
  ZEROS.zetaK20 = ZEROS.zeta.concat(ZEROS.chi20).sort(bySize);
  ZEROS.zetaK23 = ZEROS.zeta.concat(ZEROS.chi23).sort(bySize);

  /* kind: which evaluator and which gamma factor.
   * centre: the critical line. eps: the root number.
   * view: the default window. maxIm: how far the method stays trustworthy.
   * cost: relative price of one point, used to pick a grid resolution. */
  var FAMILY = {
    zeta:    { kind: 'zeta', centre: 0.5, eps: 1, view: [-2, 3, -1, 35],
               maxIm: 80, cost: 1 },
    chi3:    { kind: 'dirichlet', ch: 'chi3', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 70, cost: 2 },
    chi20:   { kind: 'dirichlet', ch: 'chi20', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 45, cost: 6 },
    chi23:   { kind: 'dirichlet', ch: 'chi23', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 45, cost: 16 },
    zetaK3:  { kind: 'dedekind', ch: 'chi3', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 70, cost: 3 },
    zetaK20: { kind: 'dedekind', ch: 'chi20', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 45, cost: 7 },
    zetaK23: { kind: 'dedekind', ch: 'chi23', centre: 0.5, eps: 1,
               view: [-2, 3, -1, 35], maxIm: 45, cost: 17 },
    /* The smoothed sum loses roughly e^{pi t / 2} worth of significance to
     * cancellation, so past |Im s| ~ 20 double precision runs out. */
    '11a1':  { kind: 'elliptic', centre: 1, eps: 1, view: [-0.5, 2.5, -2, 18],
               maxIm: 20, cost: 900 },
    '37a1':  { kind: 'elliptic', centre: 1, eps: -1, view: [-0.5, 2.5, -2, 18],
               maxIm: 20, cost: 1600 }
  };

  function evaluate(id, sr, si) {
    var f = FAMILY[id];
    switch (f.kind) {
      case 'zeta':      return zeta(sr, si);
      case 'dirichlet': return dirichletL(sr, si, CHARS[f.ch]);
      case 'dedekind':  return dedekind(sr, si, CHARS[f.ch]);
      default:          return ellipticL(sr, si, CURVES[id]);
    }
  }

  /* ------------------------------------------------------- Z-function --- */
  /* The completed function Lambda = G * L satisfies Lambda(s) = eps
   * Lambda(w-s), which forces Lambda(centre + it)/sqrt(eps) to be real. Divide
   * by |G| and you get a real function of t whose modulus is exactly
   * |L(centre+it)| and whose sign changes are exactly the zeros on the line.
   * This is Hardy's Z for zeta, and it is what LMFDB plots. */

  /* the gamma factor G(s), including the conductor power */
  function gammaFactor(id, sr, si) {
    var f = FAMILY[id], g, h;
    if (f.kind === 'zeta') {
      g = rpow(Math.PI, -sr / 2, -si / 2);
      h = cgamma(sr / 2, si / 2);
      return cmul(g[0], g[1], h[0], h[1]);
    }
    if (f.kind === 'dirichlet' || f.kind === 'dedekind') {
      /* all our characters are odd, so the shift is (s+1)/2 */
      var q = CHARS[f.ch].q;
      g = rpow(q / Math.PI, (sr + 1) / 2, si / 2);
      h = cgamma((sr + 1) / 2, si / 2);
      var chiG = cmul(g[0], g[1], h[0], h[1]);
      if (f.kind === 'dirichlet') return chiG;
      /* zeta_K = zeta * L(chi), and the gamma factors multiply to match */
      var zg = gammaFactor('zeta', sr, si);
      return cmul(zg[0], zg[1], chiG[0], chiG[1]);
    }
    var Nc = CURVES[id].N;
    g = rpow(Nc, sr / 2, si / 2);
    var tp = rpow(TWO_PI, -sr, -si);
    h = cgamma(sr, si);
    var a = cmul(g[0], g[1], tp[0], tp[1]);
    return cmul(a[0], a[1], h[0], h[1]);
  }

  function zfunction(id, t) {
    var f = FAMILY[id], c = f.centre;
    var G = gammaFactor(id, c, t);
    var L = evaluate(id, c, t);
    var P = cmul(G[0], G[1], L[0], L[1]);           /* Lambda(c+it) */
    /* divide by sqrt(eps): 1 when eps = +1, i when eps = -1 */
    var zr = f.eps > 0 ? P[0] : P[1];
    var mag = Math.hypot(G[0], G[1]);
    return mag === 0 ? 0 : zr / mag;
  }

  /* Z sampled on [t0, t1]; also returns |L| = |Z|, so the caller can scale. */
  function zline(id, t0, t1, n) {
    var out = new Float64Array(n);
    for (var i = 0; i < n; i++) out[i] = zfunction(id, t0 + i * (t1 - t0) / (n - 1));
    return out;
  }

  /* ------------------------------------------------------ grid filling --- */
  /* n^{-s} factors as exp(-sigma log n) -- which depends only on the column --
   * times cos/sin(t log n) -- which depend only on the row. Precomputing both
   * turns the inner loop into multiply-adds, and that is what makes zooming
   * feel immediate. The elliptic-curve sum has no such factorisation and is
   * evaluated point by point instead. */

  function powerGrid(M, x0, dx, W) {
    var pw = new Float64Array(W * M), logs = new Float64Array(M + 1), n, j;
    for (n = 1; n <= M; n++) logs[n] = Math.log(n);
    for (j = 0; j < W; j++) {
      var sg = x0 + j * dx;
      for (n = 1; n <= M; n++) pw[j * M + n - 1] = Math.exp(-sg * logs[n]);
    }
    return { pw: pw, logs: logs };
  }

  /* One component (zeta itself, or one Dirichlet L) over rows [r0, r1),
   * written into `dst` as interleaved (re, im). `ch` null means zeta. */
  function hurwitzBand(comp, x0, dx, W, ys, r0, r1, dst) {
    var q = comp.q, v = comp.v, M = comp.M, top = comp.top;
    var pw = comp.P.pw, logs = comp.P.logs;
    var iArr = comp.iArr, cArr = comp.cArr, K = iArr.length;
    var As = comp.As, Ac = comp.Ac, R = As.length;
    var Qs = comp.Qs, Qc = comp.Qc;
    var cs = comp.cs, sn = comp.sn;
    var i, j, n, a, k;

    for (i = r0; i < r1; i++) {
      var t = ys[i];
      for (n = 1; n <= top; n++) {
        var th = t * logs[n];
        cs[n] = Math.cos(th); sn[n] = Math.sin(th);
      }
      for (j = 0; j < W; j++) {
        var base = j * top, re = 0, im = 0;
        for (k = 0; k < K; k++) {
          var nn = iArr[k], p = cArr[k] * pw[base + nn];
          re += p * cs[nn + 1];
          im -= p * sn[nn + 1];
        }
        var sg = x0 + j * dx;
        var pc = poleCombo(sg, t, Qs, Qc);
        re += pc[0]; im += pc[1];

        for (a = 0; a < R; a++) {
          var Q = As[a], ca = Ac[a];
          var e = pw[base + Q - 1];
          var wr = e * cs[Q], wi = -e * sn[Q];
          re += ca * 0.5 * wr;
          im += ca * 0.5 * wi;
          var pr = sg, pi = t, f = q / Q, fp = f;
          for (var jj = 1; jj <= EM_C.length; jj++) {
            var cc = ca * EM_C[jj - 1] * fp;
            re += cc * (pr * wr - pi * wi);
            im += cc * (pr * wi + pi * wr);
            var m0 = pr * (sg + 2 * jj - 1) - pi * t;
            var m1 = pr * t + pi * (sg + 2 * jj - 1);
            pr = m0 * (sg + 2 * jj) - m1 * t;
            pi = m0 * t + m1 * (sg + 2 * jj);
            fp *= f * f;
          }
        }
        var o = 2 * (i * W + j);
        dst[o] = re; dst[o + 1] = im;
      }
    }
  }

  function makeComponent(ch, N, x0, dx, W) {
    var q = ch ? ch.q : 1;
    var v = ch ? ch.v : [1];
    var M = q * N, top = M + q, n, a;
    var idx = [], coef = [];
    for (n = 1; n <= M; n++) {
      var c = v[(n - 1) % q];
      if (c) { idx.push(n - 1); coef.push(c); }
    }
    var As = [], Ac = [], Qs = [], Qc = [];
    for (a = 1; a <= q; a++) {
      if (!v[a - 1]) continue;
      As.push(M + a); Ac.push(v[a - 1]);
      Qs.push(M + a); Qc.push(v[a - 1] / q);
    }
    return {
      q: q, v: v, M: M, top: top,
      P: powerGrid(top, x0, dx, W),
      iArr: new Int32Array(idx), cArr: new Float64Array(coef),
      As: As, Ac: Ac, Qs: Qs, Qc: Qc,
      cs: new Float64Array(top + 1), sn: new Float64Array(top + 1)
    };
  }

  /* A grid computation that can be advanced a band at a time, so a slow one
   * paints as it goes instead of freezing the page. Returns { out, H, step }. */
  function gridJob(id, x0, x1, y0, y1, W, H) {
    var f = FAMILY[id];
    var out = new Float64Array(2 * W * H);
    var dx = (x1 - x0) / (W - 1), i;
    var ys = new Float64Array(H);
    for (i = 0; i < H; i++) ys[i] = y1 - i * (y1 - y0) / (H - 1);

    if (f.kind === 'elliptic') {
      var curve = CURVES[id];
      return { out: out, H: H, step: function (r0, r1) {
        for (var a = r0; a < r1; a++) {
          for (var j = 0; j < W; j++) {
            var v = ellipticL(x0 + j * dx, ys[a], curve);
            var o = 2 * (a * W + j);
            out[o] = v[0]; out[o + 1] = v[1];
          }
        }
      } };
    }

    var maxT = Math.max(Math.abs(y0), Math.abs(y1));
    var N = Math.max(16, Math.ceil(maxT) + 16);

    if (f.kind === 'dedekind') {
      /* zeta_K = zeta * L(s, chi_D): compute both, multiply in place */
      var cz = makeComponent(null, N, x0, dx, W);
      var cc = makeComponent(CHARS[f.ch], N, x0, dx, W);
      var tmp = new Float64Array(2 * W * H);
      return { out: out, H: H, step: function (r0, r1) {
        hurwitzBand(cz, x0, dx, W, ys, r0, r1, tmp);
        hurwitzBand(cc, x0, dx, W, ys, r0, r1, out);
        for (var k = r0 * W; k < r1 * W; k++) {
          var o = 2 * k;
          var ar = tmp[o], ai = tmp[o + 1], br = out[o], bi = out[o + 1];
          out[o] = ar * br - ai * bi;
          out[o + 1] = ar * bi + ai * br;
        }
      } };
    }

    var comp = makeComponent(f.kind === 'zeta' ? null : CHARS[f.ch], N, x0, dx, W);
    return { out: out, H: H, step: function (r0, r1) {
      hurwitzBand(comp, x0, dx, W, ys, r0, r1, out);
    } };
  }

  /* Convenience: run the whole job at once. */
  function grid(id, x0, x1, y0, y1, W, H) {
    var job = gridJob(id, x0, x1, y0, y1, W, H);
    job.step(0, H);
    return job.out;
  }

  global.LF = {
    cgamma: cgamma, incGammaUpper: incGammaUpper, hurwitzZeta: hurwitzZeta,
    zeta: zeta, dirichletL: dirichletL, dedekind: dedekind,
    ellipticL: ellipticL, evaluate: evaluate,
    gammaFactor: gammaFactor, zfunction: zfunction, zline: zline,
    grid: grid, gridJob: gridJob,
    CHARS: CHARS, CURVES: CURVES, ZEROS: ZEROS, FAMILY: FAMILY
  };
})(typeof window !== 'undefined' ? window : this);
