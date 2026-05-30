#!/bin/bash

# Define log file
LOG_FILE="$(pwd)/run_log.txt"

echo "Starting run at $(date)" >> "$LOG_FILE"

# Loop through specific 30m directories
# "30m-mlt++" "30m-normal" "30m-ours" "30m-supereduction_a=2" "30m-supereduction_a=5"
# for dir in "30m-ours" "30m-mlt++" "30m-normal" "30m-supereduction_a=2" "30m-supereduction_a=5" "40m-ours" "40m-normal" "40m-mlt++" "40m-supereduction_a=2" "50m-ours" "50m-mlt++" "50m-normal" "50m-supereduction_a=2"; do
# for dir in "120m-ours" "120m-supereduction_a=2" "160m-ours" "160m-supereduction_a=2"; do
# for dir in "60m-ours_reduction" "60m-supereduction_a=2_reduction" "80m-ours_reduction" "80m-supereduction_a=2_reduction" ; do
for dir in {10,30,40,60}m-*; do
    [[ "$dir" == *wours_post || "$dir" == *_reduction || "$dir" == *" _reduction" ]] && continue
    if [ -d "$dir" ]; then
        echo "Entering $dir" >> "$LOG_FILE"
        cd "$dir" || continue

        # Run ./mk
        CMD="./mk"
        START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        echo "Running $CMD..." >> "$LOG_FILE"
        $CMD
        MK_EXIT_CODE=$?
        END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        # Log ./mk
        echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, ExitCode: $MK_EXIT_CODE" >> "$LOG_FILE"

        if [ $MK_EXIT_CODE -eq 0 ]; then
            CMD="./rn"
            START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            echo "Running $CMD..." >> "$LOG_FILE"
            $CMD
            RN_EXIT_CODE=$?
            END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            # Log ./rn
            echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, ExitCode: $RN_EXIT_CODE" >> "$LOG_FILE"
        else
            echo "Skipping ./rn because ./mk failed in $dir" >> "$LOG_FILE"
            echo "Folder: $dir, Command: ./rn, Start: -, End: -, Status: SKIPPED (mk failed)" >> "$LOG_FILE"
        fi

        # Go back to parent directory
        cd ..
    fi
done
