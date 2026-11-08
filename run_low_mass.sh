#!/bin/bash
# Usage: ./run_low_mass.sh [folder ...]
#
# Runs the low-mass (10, 20 and 30 Msun) *_test folders, which evolve to Si burning
# with the approx21_cr60_plus_co56 net, Ledoux + semiconvection and overshoot.
# Same flow as run.sh: clear photos, ./mk, then ../watchdog.sh (./rn / ./re).
#
# With no arguments all folders in DIRS below are run, one after the other.
# Example (single folder): ./run_low_mass.sh 10m-ours_reduction_post_test
#
# Progress is written to run_log_low_mass.txt (kept apart from run_log.txt).

cd "$(dirname "$0")" || exit 1

LOG_FILE="$(pwd)/run_log_low_mass.txt"

DIRS=(
    10m-ours_reduction_test
    10m-ours_reduction_post_test
    10m-supereduction_a=2_reduction_test
    10m-supereduction_a=2_reduction_post_test
    20m-ours_reduction_test
    20m-ours_reduction_post_test
    20m-supereduction_a=2_reduction_test
    20m-supereduction_a=2_reduction_post_test
    30m-ours_reduction_test
    30m-ours_reduction_post_test
    30m-supereduction_a=2_reduction_test
    30m-supereduction_a=2_reduction_post_test
    30m-normal_reduction_test
    30m-normal_reduction_post_test
    30m-mlt++_reduction_test
    30m-mlt++_reduction_post_test
)

if [ $# -gt 0 ]; then
    DIRS=("$@")
fi

echo "Starting low mass run at $(date)" >> "$LOG_FILE"

for dir in "${DIRS[@]}"; do
    if [ -d "$dir" ]; then
        echo "Entering $dir" >> "$LOG_FILE"
        cd "$dir" || continue

        mkdir -p photos LOGS
        rm -rf photos/*
        rm -f restart_photo watchdog.log
        echo "Cleared photos in $dir" >> "$LOG_FILE"

        CMD="./mk"
        START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        START_EPOCH=$(date +%s)
        echo "Running $dir $CMD..." >> "$LOG_FILE"
        $CMD
        MK_EXIT_CODE=$?
        END_EPOCH=$(date +%s)
        END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        TOTAL_SEC=$((END_EPOCH - START_EPOCH))
        TOTAL_TIME=$(printf "%02d:%02d:%02d" $((TOTAL_SEC/3600)) $((TOTAL_SEC%3600/60)) $((TOTAL_SEC%60)))
        echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, TotalTime: $TOTAL_TIME, ExitCode: $MK_EXIT_CODE" >> "$LOG_FILE"

        if [ $MK_EXIT_CODE -eq 0 ]; then
            CMD="../watchdog.sh"
            START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            START_EPOCH=$(date +%s)
            echo "Running $dir $CMD..." >> "$LOG_FILE"
            $CMD 2>&1 | tee /dev/null
            RN_EXIT_CODE=${PIPESTATUS[0]}
            END_EPOCH=$(date +%s)
            END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            TOTAL_SEC=$((END_EPOCH - START_EPOCH))
            TOTAL_TIME=$(printf "%02d:%02d:%02d" $((TOTAL_SEC/3600)) $((TOTAL_SEC%3600/60)) $((TOTAL_SEC%60)))
            TERM_CODE=$(grep -o 'termination code:.*' "watchdog.log" 2>/dev/null | tail -1 | sed 's/termination code: *//')
            echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, TotalTime: $TOTAL_TIME, ExitCode: $RN_EXIT_CODE, TermCode: $TERM_CODE" >> "$LOG_FILE"
            [ $TOTAL_SEC -lt 1200 ] && echo "===WARNING===" >> "$LOG_FILE"
        else
            echo "Skipping ../watchdog.sh because ./mk failed in $dir" >> "$LOG_FILE"
            echo "Folder: $dir, Command: ../watchdog.sh, Start: -, End: -, Status: SKIPPED (mk failed)" >> "$LOG_FILE"
        fi

        cd ..
    else
        echo "Folder not found, skipping: $dir" >> "$LOG_FILE"
    fi
done

echo "Finished low mass run at $(date)" >> "$LOG_FILE"
