#!/usr/bin/env python3
"""Check that every ordinate tabulated in the app is reproduced by the Magma runs.

The app's zero tables are taken from Odlyzko and LMFDB; the Magma scripts here
re-derive them. It is easy for those two to drift apart silently -- a table
extended without extending the corresponding scan height, say -- so this
compares them directly and fails loudly.

    python3 scripts/verify-tables.py

It reads the zero lists straight out of app/lfunctions.js and the ZERO lines
out of results/*.txt, so it checks what is actually shipped rather than a
separate copy of the numbers.

Note what this does and does not establish. Agreement here means the tabulated
ordinates are reproduced by an independent computation. It is not a proof that
either list is complete: both come from sign-change searches on the critical
line, which cannot see a zero of even order, two zeros inside one step, or any
zero off the line.
"""
import os
import re
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
JS = os.path.join(ROOT, "app", "lfunctions.js")
RESULTS = os.path.join(ROOT, "results")

# app key -> (Magma tag, results file)
PAIRS = [
    ("zeta",  "ZETA",  "lfunction-data-20260921.txt"),
    ("chi3",  "CHI3",  "lfunction-data-20260921.txt"),
    ("chi20", "D20",   "quad-char-zeros-20260921.txt"),
    ("chi23", "D23",   "quad-char-zeros-20260921.txt"),
    ("11a1",  "E11A1", "lfunction-data-20260921.txt"),
    ("37a1",  "E37A1", "lfunction-data-20260921.txt"),
]

TOL = 1e-9


def tabulated(js, key):
    m = re.search(re.escape(key) + r"'?\s*:\s*\[([^\]]*)\]", js)
    if not m:
        sys.exit("no zero list for %r in app/lfunctions.js" % key)
    return sorted(float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip())


def computed(path, tag):
    out = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(tag + " ZERO "):
                out.append(float(line.split()[2]))
    return sorted(out)


def main():
    js = open(JS, encoding="utf-8").read()
    bad = []
    print("%-7s %9s %6s %11s %10s" % ("", "tabulated", "magma", "scanned to", "max diff"))
    for key, tag, fn in PAIRS:
        path = os.path.join(RESULTS, fn)
        if not os.path.exists(path):
            bad.append("%s: missing %s" % (key, fn))
            continue
        tab = [x for x in tabulated(js, key) if x > 0]
        mag = computed(path, tag)
        if not mag:
            bad.append("%s: no '%s ZERO' lines in %s" % (key, tag, fn))
            continue
        missing = [x for x in tab if min(abs(x - y) for y in mag) > TOL]
        worst = max((min(abs(x - y) for y in mag) for x in tab), default=0.0)
        print("%-7s %9d %6d %11.3f %10.2e %s"
              % (key, len(tab), len(mag), max(mag), worst,
                 "" if not missing else "<-- NOT REPRODUCED"))
        for x in missing:
            bad.append("%s: %.9f is tabulated but not found by Magma "
                       "(scan reaches %.3f)" % (key, x, max(mag)))

    print()
    if bad:
        print("FAILED")
        for b in bad:
            print("  " + b)
        sys.exit(1)
    print("every tabulated ordinate is reproduced by the Magma runs")


if __name__ == "__main__":
    main()
