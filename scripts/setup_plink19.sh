#!/usr/bin/env bash
set -euo pipefail

PLINK_VERSION="${PLINK_VERSION:-20230116}"
PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip"
PLINK_ROOT="${PLINK_ROOT:-$HOME/.local/plink/1.9}"
BIN_DIR="${BIN_DIR:-$HOME/.local/bin}"

mkdir -p "$PLINK_ROOT" "$BIN_DIR"

tmp_dir="$(mktemp -d)"
cleanup() {
  rm -rf "$tmp_dir"
}
trap cleanup EXIT

echo "Downloading PLINK 1.9 (${PLINK_VERSION})..."
if command -v curl >/dev/null 2>&1; then
  curl -L -o "$tmp_dir/plink.zip" "$PLINK_URL"
elif command -v wget >/dev/null 2>&1; then
  wget -O "$tmp_dir/plink.zip" "$PLINK_URL"
else
  echo "Error: curl or wget is required to download PLINK." >&2
  exit 1
fi

unzip -o "$tmp_dir/plink.zip" -d "$tmp_dir" >/dev/null

plink_path="$tmp_dir/plink"
if [[ ! -f "$plink_path" ]]; then
  plink_path="$(find "$tmp_dir" -maxdepth 2 -type f -name plink | head -n 1 || true)"
fi

if [[ -z "${plink_path:-}" || ! -f "$plink_path" ]]; then
  echo "Error: could not locate plink binary after unzip." >&2
  exit 1
fi

mv "$plink_path" "$PLINK_ROOT/plink"
chmod +x "$PLINK_ROOT/plink"
ln -sf "$PLINK_ROOT/plink" "$BIN_DIR/plink"

echo "PLINK installed at: $PLINK_ROOT/plink"

# Ensure ~/.local/bin is on PATH permanently (persists across Cloud Shell sessions)
if ! grep -q "$BIN_DIR" "$HOME/.bashrc" 2>/dev/null; then
  echo "export PATH=\"$BIN_DIR:\$PATH\"" >> "$HOME/.bashrc"
  echo "Added $BIN_DIR to ~/.bashrc (will persist across sessions)."
fi
export PATH="$BIN_DIR:$PATH"
echo "plink is ready. Run 'plink --help' to verify."
