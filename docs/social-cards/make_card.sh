#!/usr/bin/env bash
# ── Build a 1200x630 social card: book cover (left) + case-study figure (right) ──
#
# Usage: social-cards/make_card.sh <dir[/page]> <figure> [docs_dir]
#
#   <dir[/page]>  case-study directory, e.g. roaches. If the page name
#                 differs from the directory (dogs has two pages), give
#                 it as dir/page, e.g. dogs/dogs_stan.
#   <figure>      either a figure NUMBER (N-th numbered figure in the
#                 rendered page, in document order) or a fig- chunk
#                 label (with or without the trailing "-1.png").
#   [docs_dir]    where the rendered output lives (default:
#                 ../BW-gh-pages/docs relative to this repo)
#
# Output: social-cards/<page>.png  (committed in the main branch;
#         referenced from the case study front matter via image:)
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
COVER="$REPO_DIR/9780367490188.png"          # full cover (narrower than the crop)
OUT_DIR="$REPO_DIR/social-cards"
DOCS_DIR="${3:-$REPO_DIR/../BW-gh-pages/docs}"

arg="$1"; fig="$2"
dir="${arg%%/*}"; page="${arg##*/}"          # split dir/page; page==dir if no slash

html="$DOCS_DIR/$dir/$page.html"
figdir="$DOCS_DIR/$dir/${page}_files/figure-html"

mkdir -p "$OUT_DIR"
out="$OUT_DIR/$page.png"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# ── Resolve the figure into $tmp/figraw.png ──
# embed-resources inlines every figure as base64, so the most uniform
# approach is to pull the N-th base64 PNG (document order) from the HTML.
# A non-numeric <figure> is treated as a fig- file name in the _files dir.
if [[ "$fig" =~ ^[0-9]+$ ]]; then
  [ -f "$html" ] || { echo "No rendered HTML: $html" >&2; exit 1; }
  grep -o 'data:image/png;base64,[A-Za-z0-9+/=]\+' "$html" \
    | sed -n "${fig}p" | sed 's/^data:image\/png;base64,//' \
    | base64 -d > "$tmp/figraw.png" 2>/dev/null || true
  [ -s "$tmp/figraw.png" ] || { echo "Page $page has fewer than $fig figures" >&2; exit 1; }
  figdesc="figure #$fig"
else
  if   [ -f "$figdir/$fig" ];         then figpath="$figdir/$fig"
  elif [ -f "$figdir/${fig}.png" ];   then figpath="$figdir/${fig}.png"
  elif [ -f "$figdir/${fig}-1.png" ]; then figpath="$figdir/${fig}-1.png"
  else
    echo "Figure '$fig' not found in $figdir" >&2
    echo "Available:" >&2; ls "$figdir"/*.png 2>/dev/null | xargs -n1 basename >&2
    exit 1
  fi
  cp "$figpath" "$tmp/figraw.png"
  figdesc="$(basename "$figpath")"
fi

# ── Left panel: full cover scaled to 560px height (~409px wide) ──
convert "$COVER" -resize x560 "$tmp/cover.png"

# ── Right panel: figure fit into a 680x560 box on white, hairline border ──
convert "$tmp/figraw.png" -resize 680x560 -background white -gravity center \
  -extent 680x560 -bordercolor "#dddddd" -border 1 "$tmp/fig.png"

# ── Compose on a 1200x630 white canvas: cover left, figure right ──
convert -size 1200x630 xc:white \
  "$tmp/cover.png" -gravity west -geometry +50+0 -composite \
  "$tmp/fig.png"   -gravity east -geometry +40+0 -composite \
  -strip "$out"

echo "Wrote $out  ($figdesc)"
