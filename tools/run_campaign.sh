#!/usr/bin/env bash
# Re-run the complete openEMS verification campaign of the technical report.
#
#   bash tools/run_campaign.sh [campaign_dir]        (default: campaign/, git-ignored)
#
# Then rebuild everything from the new results:
#   python docs/report/collect.py campaign           # results/openems/{studies,reference,...}
#   python tools/verify_results.py campaign          # independent check of every number
#   python docs/report/build_report.py               # the PDF
#   python docs/report/update_readme.py              # README results block
#
# 24 runs. On two CPU cores the whole campaign takes about 2 hours; the
# h = 0.0125 mm reference run alone takes about 25 minutes and needs ~2 GB RAM.
# Linux or WSL with openEMS installed (see openems/README.md).
set -u
cd "$(dirname "$0")/.."
C="${1:-campaign}"
mkdir -p "$C"
LOG="$C/campaign.log"

run() {                                   # run <name> [simulate_openems.py options...]
  local name=$1; shift
  if [ -f "$C/$name/run_info.json" ]; then echo "SKIP $name (done)" | tee -a "$LOG"; return; fi
  echo "START $name $(date +%T)" | tee -a "$LOG"
  if python openems/simulate_openems.py --out "$C/$name" --sim-dir "$C/_sim/$name" "$@" > "$C/$name.log" 2>&1; then
    echo "DONE  $name  $(grep -E 'coupling  |directivity  ' "$C/$name.log" | tr -s ' ' | tr '\n' ' ')" | tee -a "$LOG"
  else
    echo "FAIL  $name (see $C/$name.log)" | tee -a "$LOG"
  fi
  rm -rf "$C/_sim/$name"
}
variant() {                               # variant <name> [moreno.design options...]
  local name=$1; shift
  python -m moreno.design "$@" --out "$C/$name.json" > /dev/null
}

# design variants for the chart and the tolerance study
for L in 1.6 1.8 2.2 2.4 2.6; do variant "L$L" --L "$L"; done
variant t020 --t 0.20; variant t030 --t 0.30
variant w030 --w 0.30; variant w040 --w 0.40

# 1. slot-mesh convergence, 1/3-2/3 edge rule (Table 4)
run c050  --fine 0.05
run c035  --fine 0.035
run c025  --fine 0.025
run c0175 --fine 0.0175
run c0125 --fine 0.0125
# 2. the same with mesh lines on the slot edges: brackets the limit from the other side
run s035  --staircase --fine 0.035
run s025  --staircase --fine 0.025
run s0175 --staircase --fine 0.0175
# 3. other numerical parameters, each changed alone (Table 5)
run n_res30  --fine 0.035 --res 30
run n_feed25 --fine 0.035 --feed 2.5
run n_e4     --fine 0.035 --end-criteria 1e-4
# 4. full 4x4 scattering matrix (Table 3) and field maps (Fig. 6)
run full   --fine 0.035 --excite all
run fields --fine 0.035 --dump-fields
# 5. design chart (Table 7) and its mesh-shift check at h = 0.025 mm
for L in 1.6 1.8 2.2 2.4 2.6; do run "L$L" --fine 0.035 --design "$C/L$L.json" --max-steps 150000; done
run o025_L1.6 --fine 0.025 --design "$C/L1.6.json"
run o025_L2.4 --fine 0.025 --design "$C/L2.4.json" --max-steps 150000
# 6. manufacturing tolerances (Table 8)
for v in t020 t030 w030 w040; do run "$v" --fine 0.035 --design "$C/$v.json"; done

echo "CAMPAIGN DONE $(date +%T)" | tee -a "$LOG"
