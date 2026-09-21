# Data for rh-exploration-app

This directory is gitignored. Record here everything needed to find or
regenerate the real files.

| File | Size | Where it actually lives | How it was produced | sha256 |
|---|---|---|---|---|
| `odlyzko-zeros2.txt` | 105,979 B | <https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros2> | Downloaded 2026-09-21. Andrew Odlyzko's table of the first 100 zeros of ζ(s), each to over 1000 decimal places. | `0439d90a4c025d1ab3ed25f2241f27afeb6d01e651d95672267783b859ee170f` |

## Format of `odlyzko-zeros2.txt`

One zero per record. A record starts on a line beginning with two spaces and
the integer part; continuation lines carry 66 more digits each and begin with a
single space. To get one number per line:

```bash
awk '
  /^  [0-9]/ { if (buf != "") print buf; buf="" }
  /^ *[0-9]/ { s=$0; gsub(/ /,"",s); buf = buf s }
  END { if (buf != "") print buf }
' data/odlyzko-zeros2.txt
```

The first 20 of these, truncated to double precision, are what
`app/lfunctions.js` carries in `ZEROS.zeta`. They agree to every digit with the
values computed independently in `scripts/lfunction-data.m`.

## Zeros not from this file

The other zero lists in `app/lfunctions.js` are not large enough to need a data
file, and are recorded in `results/` instead:

- **L(s, χ₋₃)** — from `scripts/lfunction-data.m`, agreeing with
  [LMFDB 1-3-3.2-r1-0-0](https://www.lmfdb.org/L/Character/Dirichlet/3/2/) to
  the three decimals LMFDB prints.
- **L(E, s) for 11a1 and 37a1** — from
  [LMFDB 2-11-1.1-c1-0-0](https://www.lmfdb.org/L/EllipticCurve/Q/11/a/) and
  [2-37-1.1-c1-0-1](https://www.lmfdb.org/L/EllipticCurve/Q/37/a/), which give
  about 30 digits, cross-checked against `scripts/lfunction-data.m`.
