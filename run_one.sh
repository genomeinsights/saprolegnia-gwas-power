#!/usr/bin/env bash
## run_one.sh C REP -- SLiM -> recapitate/mutate -> two-stage GWAS for one replicate.
## Results: out/res/C<C>_r<REP>.res.rds ; logs: out/logs/. Large files are deleted.
set -euo pipefail
cd "$(dirname "$0")"
C="$1"; REP="$2"
PY="${PY:-python3}"   # needs msprime >= 1.3 and tskit
ID="C${C}_r${REP}"
SEED=$(( REP * 1000 + $(printf '%.0f' "$(echo "$C * 10000" | bc)") % 1000 ))
mkdir -p out/tmp out/res out/logs
DEPTH_HIGH="${DEPTH_HIGH:-15}"; export DEPTH_HIGH
RESDIR=out/res; [ "$DEPTH_HIGH" != 15 ] && RESDIR="out/res_d${DEPTH_HIGH}"
mkdir -p "$RESDIR"
RES="${RESDIR}/${ID}.res.rds"
[ -f "$RES" ] && exit 0
TMP="out/tmp/${ID}_d${DEPTH_HIGH}"
{
  slim -d C="$C" -d SEED="$SEED" -d "OUT='${TMP}'" slim/clonal_metapop.slim | tail -1
  "$PY" -W ignore py/ts_to_geno.py --trees "${TMP}.trees" --samples "${TMP}.samples.csv" \
        --out "$TMP" --seed "$SEED"
  Rscript R/analyse_rep.R "$TMP" "$SEED" "$RES"
} > "out/logs/${ID}_d${DEPTH_HIGH}.log" 2>&1 || { echo "FAILED $ID"; rm -f "${TMP}".*; exit 1; }
rm -f "${TMP}".*
echo "ok $ID"
