# The Riemann Hypothesis — an explorer

**[Open the app →](https://alozanoroble.github.io/riemann-hypothesis-explorer/)**

A teaching tool for the Riemann Hypothesis and its three standard
generalisations. Each one gets the same treatment: the conjecture stated
precisely, the *L*-function defined, its coefficients and Euler product, the
critical strip and functional equation, a table of zeros, a live plot, and
references.

| | Conjecture | Example *L*-functions |
|---|---|---|
| **RH** | zeros of ζ(*s*) | the Riemann zeta function |
| **ERH** | zeros of ζ<sub>*K*</sub>(*s*) | ℚ(√−5) with *h* = 2, ℚ(√−23) with *h* = 3, ℚ(√−3) with *h* = 1 |
| **GRH** | zeros of *L*(*s*, χ) | quadratic characters of conductor 3, 20, 23 |
| **Grand RH** | zeros of automorphic *L* | elliptic curves 11a1 (rank 0) and 37a1 (rank 1) |

Every plot has two views. The **complex plane** view colours the strip by
|*L*(*s*)|, white at the zeros, so they appear as a row of white points on a
single vertical line. The **critical line** view plots the *Z*-function — the
completed *L*-function rotated so that it is real — where the zeros become
crossings, the way LMFDB draws them.

Everything is computed live in the browser. Zoom with the scroll wheel, drag to
pan, hover to read off a value.

## Where the numbers come from

Nothing here is quoted from memory. Every zero comes from a published table or
from a computation in this repository, and usually both:

- **ζ(*s*)** — [Odlyzko's table](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/)
  of the first 100 zeros to over 1000 decimal places.
- **Dirichlet *L*-functions** — LMFDB
  [1-3-3.2-r1-0-0](https://www.lmfdb.org/L/Character/Dirichlet/3/2/),
  [1-20-20.19-r1-0-0](https://www.lmfdb.org/L/Character/Dirichlet/20/19/),
  [1-23-23.22-r1-0-0](https://www.lmfdb.org/L/Character/Dirichlet/23/22/).
- **Elliptic curves** — LMFDB
  [2-11-1.1-c1-0-0](https://www.lmfdb.org/L/EllipticCurve/Q/11/a/) and
  [2-37-1.1-c1-0-1](https://www.lmfdb.org/L/EllipticCurve/Q/37/a/).
- **Dedekind zeta** — assembled from its two factors, since
  ζ<sub>*K*</sub> = ζ · *L*(*s*, χ<sub>*D*</sub>) for a quadratic field.

The Magma scripts in `scripts/` reproduce all of them independently, and agree
with the published values to every digit those sources print. `results/` holds
the output they produced.

## Computing the zeros yourself

The app shows the code, and it is here too:

- [`scripts/lfunction-data.m`](scripts/lfunction-data.m) — Magma: zeros,
  coefficients and reference values for ζ, a Dirichlet *L*-function and two
  elliptic curves.
- [`scripts/dedekind-data.m`](scripts/dedekind-data.m),
  [`scripts/quad-char-values.m`](scripts/quad-char-values.m),
  [`scripts/quad-char-zeros.m`](scripts/quad-char-zeros.m) — the Dedekind zeta
  material.
- [`scripts/lfunction-zeros.sage`](scripts/lfunction-zeros.sage) — the Sage
  equivalent, using Rubinstein's `lcalc`.

Magma 2.29 has no zero-finding intrinsic for the `LSer` type — `Zeros` only
accepts curves, function fields and number fields — so these scripts locate
zeros themselves. If *L* is self-dual with root number ε and motivic weight
*m*, then Λ(*s*) = ε Λ(*m*+1−*s*) forces Λ(*c* + *it*)/√ε to be real on the
critical line, and its sign changes are the zeros. A scan plus the Illinois
root-finder does the rest.

One trap worth knowing: `Sign(L)` returns 0 until you have called
`CheckFunctionalEquation(L)`, which then divides by zero further down.

> **Note on the Sage script.** It follows the documented examples in the Sage
> reference manual but has not been executed — the machine this was built on
> has Magma and no Sage. The Magma scripts are the ones that produced the
> numbers here. Corrections welcome.

## Numerical honesty

Two limits are stated in the app rather than hidden, and they are worth
knowing before you trust a picture:

- **The elliptic-curve evaluator loses precision as you go up.** It uses the
  incomplete-gamma smoothed sum, which suffers cancellation of roughly
  e<sup>π*t*/2</sup>. At Im(*s*) = 20 only about five significant digits
  survive, which is why that window stops at 18. The zeros marked on those
  plots come from LMFDB, not from this code. The other evaluators use
  Euler–Maclaurin on the Hurwitz zeta function and hold full double precision
  throughout their windows.
- **The plot axes are not to equal scale.** A 36-unit stretch of a strip five
  units wide cannot be drawn isotropically and still fit on a screen, so the
  vertical axis is compressed by about 5×. The app prints the factor.

`app/test-lfunctions.html` is the check: it compares the browser's arithmetic
against 30 Magma reference values, verifies the analytic class number formula
2π*h*/(*w*√|*D*|) for all three fields, confirms that every tabulated zero
really is a zero, that |*Z*| = |*L*| and that *Z* changes sign at each zero,
and that the grid renderer agrees with pointwise evaluation.

## Running it locally

`docs/index.html` is a single self-contained file — the whole app, no build
step, no dependencies beyond Google Fonts. Open it in a browser, or put it on
any web server.

To work on the sources in `app/`, serve them (the page will not load its script
from a `file://` URL):

```bash
python3 scripts/serve.py 8766          # serves ./app
python3 scripts/build-standalone.py    # regenerates docs/index.html
```

Use `scripts/serve.py` rather than `python3 -m http.server`: the stock server
sends no charset, and the page then renders as windows-1252 with every Greek
letter mangled.

Edit `app/index.html` and `app/lfunctions.js`; `docs/index.html` is generated,
so rebuild it rather than editing it.

## Layout

| Path | |
|---|---|
| `docs/index.html` | the built, self-contained app (what GitHub Pages serves) |
| `app/index.html` | the page source |
| `app/lfunctions.js` | the numerics |
| `app/test-lfunctions.html` | validation against Magma |
| `scripts/` | Magma, Sage, the build script and the dev server |
| `results/` | raw Magma output behind the tables |
| `data/README.md` | provenance of Odlyzko's table |

## Adding another example

For the Euler–Maclaurin families it is cheap. A new quadratic character needs
its conductor, its values on 1…*q*, and its zeros, added to `CHARS`, `ZEROS`
and `FAMILY` in `app/lfunctions.js`. A new elliptic curve needs its conductor,
root number and *a<sub>n</sub>*. One gap: every character currently included is
odd, so a real quadratic field would need the even branch of `gammaFactor`
written.

---

Álvaro Lozano-Robledo, Department of Mathematics, University of Connecticut.
