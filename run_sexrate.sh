#!/usr/bin/env bash
## run_sexrate.sh C N REP -- SLiM (deme size N) -> genotypes -> LD summaries (R/sexrate_rep.R).
set -euo pipefail
cd "$(dirname "$0")"
C="$1"; N="$2"; REP="$3"
PY="${PY:-python3}"
ID="C${C}_N${N}_r${REP}"
SEED=$(( REP * 1000 + $(printf '%.0f' "$(echo "$C * 10000" | bc)") % 1000 + N / 100 ))
mkdir -p out/sexrate/tmp out/sexrate/res out/sexrate/logs
RES="out/sexrate/res/${ID}.rds"; [ -f "$RES" ] && exit 0
TMP="out/sexrate/tmp/${ID}"
{
  slim -d C="$C" -d N="$N" -d SEED="$SEED" -d "OUT='${TMP}'" slim/clonal_metapop.slim | tail -1
  "$PY" -W ignore py/ts_to_geno.py --trees "${TMP}.trees" --samples "${TMP}.samples.csv" --out "$TMP" --seed "$SEED" --n "$N"
  Rscript R/sexrate_rep.R "$TMP" "$SEED" "$RES"
} > "out/sexrate/logs/${ID}.log" 2>&1 || { echo "FAILED $ID"; rm -f "${TMP}".*; exit 1; }
rm -f "${TMP}".*
echo "ok $ID"
