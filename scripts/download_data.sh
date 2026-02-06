#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/download_data.sh <URL> [output_file]

Notes:
- For Dropbox links, the script converts dl=0 to dl=1 for direct download.
- If output_file is omitted, the filename is inferred from the URL.
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

url="$1"
out="${2:-}"

if [[ "$url" == *"dropbox.com"* ]]; then
  if [[ "$url" == *"dl=0"* ]]; then
    url="${url//dl=0/dl=1}"
  elif [[ "$url" != *"dl=1"* ]]; then
    if [[ "$url" == *"?"* ]]; then
      url="${url}&dl=1"
    else
      url="${url}?dl=1"
    fi
  fi
fi

if [[ -z "$out" ]]; then
  out="$(basename "${url%%\?*}")"
fi

echo "Downloading to: $out"
if command -v curl >/dev/null 2>&1; then
  curl -L -o "$out" "$url"
elif command -v wget >/dev/null 2>&1; then
  wget -O "$out" "$url"
else
  echo "Error: curl or wget is required to download data." >&2
  exit 1
fi
