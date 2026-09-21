#!/usr/bin/env python3
"""Static server for local testing of the app.

Python's stock http.server sends text/html with no charset, so the browser
falls back to windows-1252 and every Greek letter in the page turns to
mojibake. This one declares UTF-8.

    python3 scripts/serve.py [port] [dir]   # dir defaults to ./app
"""
import http.server, socketserver, sys, os

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8766
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "app")


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=ROOT, **kw)

    def guess_type(self, path):
        t = super().guess_type(path)
        if t in ("text/html", "application/javascript", "text/javascript"):
            return t + "; charset=utf-8"
        return t


socketserver.TCPServer.allow_reuse_address = True
with socketserver.TCPServer(("127.0.0.1", PORT), Handler) as s:
    print("serving %s on http://127.0.0.1:%d" % (os.path.realpath(ROOT), PORT))
    s.serve_forever()
