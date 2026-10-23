#!/usr/bin/env bash
## run_loh.sh C REP [LOH_LEN] -- SLiM -> genotypes + true LOH tracts -> LOH proof of concept.
## Results: out/loh/res/C<C>_L<LOH_LEN>_r<REP>.rds ; intermediate files are deleted.
set -euo pipefail
cd "$(dirname "$0")"
C="$1"; REP="$2"; LEN="${3:-100000}"
PY="${PY:-python3}"   # needs msprime >= 1.3 and tskit
ID="C${C}_L${LEN}_r${REP}"
SEED=$(( REP * 1000 + $(printf '%.0f' "$(echo "$C * 10000" | bc)") % 1000 ))
mkdir -p out/loh/tmp out/loh/res out/loh/logs
RES="out/loh/res/${ID}.rds"
[ -f "$RES" ] && exit 0
TMP="out/loh/tmp/${ID}"
{
  slim -d C="$C" -d SEED="$SEED" -d LOH_LEN="$LEN" -d "OUT='${TMP}'" slim/clonal_metapop.slim | tail -1
  "$PY" -W ignore py/ts_to_geno.py --trees "${TMP}.trees" --samples "${TMP}.samples.csv" \
        --out "$TMP" --seed "$SEED" --loh
  Rscript R/loh_poc.R "$TMP" "$SEED" "$RES"
} > "out/loh/logs/${ID}.log" 2>&1 || { echo "FAILED $ID"; rm -f "${TMP}".*; exit 1; }
rm -f "${TMP}".*
echo "ok $ID"
