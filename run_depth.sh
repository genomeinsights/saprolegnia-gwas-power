#!/usr/bin/env bash
## run_depth.sh DEPTH C REP -- run_one.sh with the representatives sequenced at DEPTH (x).
DEPTH_HIGH="$1" exec "$(dirname "$0")/run_one.sh" "$2" "$3"
