#!/bin/bash
# Usage: ./run_one.sh <directory_name>
# Example: ./run_one.sh 10m-ours_reduction

if [ -z "$1" ]; then
    echo "Error: Please provide a directory name as parameter"
    echo "Usage: $0 <directory_name>"
    exit 1
fi

LOG_FILE="$(pwd)/run_log.txt"

echo "Starting run at $(date)" >> "$LOG_FILE"

DIRS=("$1")

for dir in "${DIRS[@]}"; do
    if [ -d "$dir" ]; then
        echo "Entering $dir" >> "$LOG_FILE"
        cd "$dir" || continue

        CMD="./mk"
        START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        START_EPOCH=$(date +%s)
        echo "Running $CMD..." >> "$LOG_FILE"
        $CMD
        MK_EXIT_CODE=$?
        END_EPOCH=$(date +%s)
        END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        TOTAL_SEC=$((END_EPOCH - START_EPOCH))
        TOTAL_TIME=$(printf "%02d:%02d:%02d" $((TOTAL_SEC/3600)) $((TOTAL_SEC%3600/60)) $((TOTAL_SEC%60)))
        echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, TotalTime: $TOTAL_TIME, ExitCode: $MK_EXIT_CODE" >> "$LOG_FILE"

        if [ $MK_EXIT_CODE -eq 0 ]; then
            CMD="./rn"
            START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            START_EPOCH=$(date +%s)
            echo "Running $CMD..." >> "$LOG_FILE"
            TMP_FILE=$(mktemp)
            $CMD 2>&1 | tee "$TMP_FILE"
            RN_EXIT_CODE=${PIPESTATUS[0]}
            END_EPOCH=$(date +%s)
            END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
            TOTAL_SEC=$((END_EPOCH - START_EPOCH))
            TOTAL_TIME=$(printf "%02d:%02d:%02d" $((TOTAL_SEC/3600)) $((TOTAL_SEC%3600/60)) $((TOTAL_SEC%60)))
            TERM_CODE=$(grep -o 'termination code:.*' "$TMP_FILE" | tail -1 | sed 's/termination code: *//')
            rm -f "$TMP_FILE"
            echo "Folder: $dir, Command: $CMD, Start: $START_TIME, End: $END_TIME, TotalTime: $TOTAL_TIME, ExitCode: $RN_EXIT_CODE, TermCode: $TERM_CODE" >> "$LOG_FILE"
            [ $TOTAL_SEC -lt 1200 ] && echo "===WARNING===" >> "$LOG_FILE"
        else
            echo "Skipping ./rn because ./mk failed in $dir" >> "$LOG_FILE"
            echo "Folder: $dir, Command: ./rn, Start: -, End: -, Status: SKIPPED (mk failed)" >> "$LOG_FILE"
        fi

        cd ..
    fi
done
