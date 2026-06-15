#!/usr/bin/env python3
# ── Add `#' image: ../social-cards/<page>.png` to each case study's front matter ──
# Reads social-cards/cards.tsv; idempotent (skips files already wired).
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

def entries():
    with open(os.path.join(HERE, "cards.tsv")) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            page = line.split("\t")[0].strip()
            yield page                       # e.g. "bioassay" or "dogs/dogs_stan"

for page in entries():
    d = page.split("/")[0]                   # directory
    p = page.split("/")[-1]                  # page / .R basename / card name
    rfile = os.path.join(REPO, d, p + ".R")
    if not os.path.exists(rfile):
        print("missing .R:", rfile); continue
    text = open(rfile).read()
    if re.search(r"^#' image:", text, re.M):
        print("skip (already wired):", page); continue
    new, n = re.subn(r"(?m)^(#' ---\n)",
                     r"\1#' image: ../social-cards/%s.png\n" % p, text, count=1)
    if n == 0:
        print("NO front-matter opener:", page); continue
    open(rfile, "w").write(new)
    print("wired:", page, "->", p + ".png")
