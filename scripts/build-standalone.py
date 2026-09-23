#!/usr/bin/env python3
"""Build a standalone, self-contained copy of the app.

app/index.html is written for the Artifact runtime, which wraps it in a
doctype, a head with charset and viewport, and a small reset, and which serves
app/lfunctions.js alongside it. Anywhere else -- GitHub Pages, a department
web server, a USB stick -- none of that happens, so this script supplies the
wrapper and inlines the script, producing one file that depends on nothing but
Google Fonts.

    python3 scripts/build-standalone.py            # -> docs/index.html

The single source of truth stays app/index.html. Do not edit docs/index.html;
rebuild it.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "app", "index.html")
JS = os.path.join(ROOT, "app", "lfunctions.js")
OUT = os.path.join(ROOT, "docs", "index.html")

# What the Artifact runtime would otherwise provide.
HEAD_OPEN = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="description" content="The Riemann Hypothesis and its generalizations \
-- ERH for Dedekind zeta functions, GRH for Dirichlet L-functions, and the Grand \
Riemann Hypothesis -- each with its L-function, its zeros, and live plots of the \
critical strip and the critical line.">
<style>
  :root {
    color-scheme: light dark;
    padding-top: env(safe-area-inset-top, 0px);
    padding-bottom: env(safe-area-inset-bottom, 0px);
  }
  img { max-width: 100%; }
  [hidden] { display: none !important; }
</style>
"""

HEAD_CLOSE = """</head>
<body>
"""

# Shared top banner linking back to alozanoroble.github.io (lives in the
# website-build repo).
TAIL = """
<script src="https://alozanoroble.github.io/banner.js" defer></script>
</body>
</html>
"""


def main():
    with open(SRC, encoding="utf-8") as f:
        page = f.read()
    with open(JS, encoding="utf-8") as f:
        js = f.read()

    # Everything up to and including the first </style> is head material
    # (title, font links, the stylesheet); the rest is the document body.
    split = page.find("</style>")
    if split < 0:
        sys.exit("no </style> found in app/index.html -- has its shape changed?")
    split += len("</style>")
    head, body = page[:split], page[split:]

    # Inline the external script.
    pattern = re.compile(r'<script src="lfunctions\.js"></script>')
    if not pattern.search(body):
        sys.exit("no <script src=\"lfunctions.js\"> found -- has its shape changed?")
    body = pattern.sub(lambda m: "<script>\n" + js + "\n</script>", body, count=1)

    # The inlined source mentions its own filename in a comment, so look for a
    # surviving *reference* rather than the bare string.
    if 'src="lfunctions.js"' in body:
        sys.exit("a script reference to lfunctions.js survived; the build would be broken")

    out = HEAD_OPEN + head + "\n" + HEAD_CLOSE + body + TAIL

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(out)

    print("wrote %s  (%.1f KB)" % (OUT, len(out.encode("utf-8")) / 1024.0))


if __name__ == "__main__":
    main()
