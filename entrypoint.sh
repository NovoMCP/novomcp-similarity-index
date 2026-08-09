#!/usr/bin/env bash
# Dispatch: `build ...` -> build_index.py, `serve ...` -> server.py.
set -euo pipefail
cmd="${1:-serve}"; shift || true
case "$cmd" in
  build) exec python /app/build_index.py "$@" ;;
  serve) exec python /app/server.py "$@" ;;
  *)     echo "unknown command: $cmd (expected 'build' or 'serve')" >&2; exit 2 ;;
esac
