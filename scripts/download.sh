#!/usr/bin/env bash
# Download a Metabolomics Workbench study (summary, data, factors, mwTab) via the REST API.
# Usage: scripts/download.sh [STUDY_ID]   (default: ST000390)
set -euo pipefail

STUDY="${1:-ST000390}"
REST="https://www.metabolomicsworkbench.org/rest"
OUT="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$OUT"

# fetch URL DEST: download to a temp file, move into place only if non-empty.
# An empty JSON list ("[]") is discarded rather than saved.
fetch() {
    local url="$1" dest="$2" tmp
    tmp="$(mktemp "$OUT/.$(basename "$dest").XXXXXX")"
    trap 'rm -f "$tmp"' EXIT
    echo "Fetching $url ..."
    curl -fsSL --retry 3 --max-time 300 "$url" -o "$tmp"
    if [[ ! -s "$tmp" ]]; then
        echo "Error: empty response for $url" >&2
        exit 1
    fi
    if [[ "$(tr -d '[:space:]' < "$tmp")" == "[]" ]]; then
        echo "  -> no records returned; not saving $(basename "$dest")"
        rm -f "$tmp" "$dest"
    else
        mv "$tmp" "$dest"
        echo "  -> $dest ($(wc -c < "$dest" | tr -d ' ') bytes)"
    fi
    trap - EXIT
}

for what in summary data factors; do
    fetch "$REST/study/study_id/$STUDY/$what" "$OUT/${STUDY}_${what}.json"
done

# mwTab is per analysis; it carries subject metadata (Stage, Age, Race) and
# the full data matrix, including unannotated features the data endpoint omits.
analyses="$(curl -fsSL --retry 3 "$REST/study/study_id/$STUDY/analysis" \
    | grep -o '"analysis_id":"AN[0-9]*"' | cut -d'"' -f4)"
if [[ -z "$analyses" ]]; then
    echo "Error: no analysis IDs found for $STUDY" >&2
    exit 1
fi
for an in $analyses; do
    fetch "$REST/study/analysis_id/$an/mwtab/txt" "$OUT/${STUDY}_${an}_mwtab.txt"
done
