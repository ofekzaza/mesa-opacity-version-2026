#!/bin/bash
# Rerun all run.sh folders (excluding 300M): recompile with ./mk, then
# resume from the last photo with ./re. Afterwards, verify every folder
# has non-empty photo files.
#
# Usage:   ./rerun_from_photo.sh            # all run.sh folders except 300M
#          ./rerun_from_photo.sh 10m-ours_reduction   # just one folder
set -u

BASE_DIR="$(pwd)"
LOG_FILE="$BASE_DIR/rerun_from_photo_log.txt"
: > "$LOG_FILE"

DIRS=(

    10m-ours_reduction
    10m-ours_reduction_post
    10m-supereduction_a=2_reduction
    10m-supereduction_a=2_reduction_post
    30m-mlt++_reduction_default
    30m-mlt++_reduction_post_default
    30m-normal_reduction_default
    30m-normal_reduction_post_default
    30m-ours_reduction_default
    30m-ours_reduction_post_default
    30m-supereduction_a=2_reduction_default
    30m-supereduction_a=2_reduction_post_default
    20m-ours_reduction
    20m-ours_reduction_post
    20m-supereduction_a=2_reduction
    20m-supereduction_a=2_reduction_post
    30m-mlt++_reduction
    30m-mlt++_reduction_post
    30m-normal_reduction
    30m-normal_reduction_post
    30m-ours_reduction
    30m-ours_reduction_post
    30m-supereduction_a=2_reduction
    30m-supereduction_a=2_reduction_post
    100m-ours_reduction
    100m-ours_reduction_post
    100m-supereduction_a=2_reduction
    100m-supereduction_a=2_reduction_post
    40m-ours_reduction
    40m-ours_reduction_post
    40m-supereduction_a=2_reduction
    40m-supereduction_a=2_reduction_post
    50m-ours_reduction
    50m-ours_reduction_post
    50m-supereduction_a=2_reduction
    50m-supereduction_a=2_reduction_post
    60m-ours_reduction
    60m-ours_reduction_post
    60m-supereduction_a=2_reduction
    60m-supereduction_a=2_reduction_post
    70m-ours_reduction
    70m-ours_reduction_post
    70m-supereduction_a=2_reduction
    70m-supereduction_a=2_reduction_post
    80m-ours_reduction
    80m-ours_reduction_post
    80m-supereduction_a=2_reduction
    80m-supereduction_a=2_reduction_post
    90m-ours_reduction
    90m-ours_reduction_post
    90m-supereduction_a=2_reduction
    90m-supereduction_a=2_reduction_post
    20m-ours_reduction_suite
)

# Optional single-folder mode, or "continue from folder F" mode (pass
# `--from F` to resume the full sequence starting at folder F).
START_AT=""
if [ $# -ge 2 ] && [ "$1" = "--from" ]; then
    START_AT="$2"
elif [ $# -gt 0 ]; then
    DIRS=("$1")
fi

photo_count() {  # $1 = folder
    find "$1/photos" -maxdepth 1 -type f -size +0c 2>/dev/null | wc -l
}

prepare_restart_photo() {  # run inside the folder
    # Pick the photo with the largest numeric name (photo filenames are model
    # numbers).  Fall back to the most recently modified file when there are
    # no numeric ones.  (ls -t is unsafe here: old photos share one mtime,
    # so it can return the *oldest* numeric photo, e.g. 1000.)
    latest=$(ls photos/ 2>/dev/null | awk '/^[0-9]+$/ {n=$1+0; if (n>max) {max=n; name=$1}} END {print name}')
    if [ -z "$latest" ]; then
        latest=$(ls -t photos/* 2>/dev/null | head -1)
    fi
    if [ -n "$latest" ]; then
        rm -f restart_photo
        cp -f "photos/$latest" restart_photo
    fi
}

echo "Starting rerun (mk + re) at $(date)" >> "$LOG_FILE"
echo "Folders: ${#DIRS[@]}" >> "$LOG_FILE"

skip_to=0
[ -n "$START_AT" ] && skip_to=1

for dir in "${DIRS[@]}"; do
    if [ "$skip_to" = 1 ]; then
        if [ "$dir" = "$START_AT" ]; then
            skip_to=0
        else
            echo "[SKIP] $dir (before --from $START_AT)" >> "$LOG_FILE"
            continue
        fi
    fi
    if [ ! -d "$BASE_DIR/$dir" ]; then
        echo "[SKIP] $dir not found" >> "$LOG_FILE"
        continue
    fi
    cd "$BASE_DIR/$dir" || continue

    np=$(photo_count "$PWD")
    if [ "$np" -eq 0 ]; then
        echo "[SKIP] $dir has no non-empty photos -> cannot ./re" >> "$LOG_FILE"
        cd "$BASE_DIR"
        continue
    fi

    # Pick the newest photo by name (model number).  This is passed to
    # ./re explicitly, because `re` chooses its own photo via `ls -tp`,
    # which tie-breaks to the OLDEST numeric photo (1000) when all photos
    # share one mtime.
    PHOTO=$(ls photos/ 2>/dev/null | awk '/^[0-9]+$/ {n=$1+0; if (n>max) {max=n; name=$1}} END {print name}')
    if [ -z "$PHOTO" ]; then
        PHOTO=$(ls -t photos/* 2>/dev/null | head -1 | xargs basename)
    fi
    if [ -z "$PHOTO" ]; then
        echo "[SKIP] $dir could not determine a photo to resume from" >> "$LOG_FILE"
        cd "$BASE_DIR"
        continue
    fi

    echo "[RUN] $dir ./mk" >> "$LOG_FILE"
    START=$(date "+%Y-%m-%d %H:%M:%S")
    ./mk > mk.out 2>&1
    MK=$?
    END=$(date "+%Y-%m-%d %H:%M:%S")
    echo "  mk start=$START end=$END exit=$MK" >> "$LOG_FILE"
    if [ "$MK" -ne 0 ]; then
        echo "  SKIP ./re because ./mk failed (see $dir/mk.out)" >> "$LOG_FILE"
        cd "$BASE_DIR"
        continue
    fi

    echo "[RUN] $dir ./re $PHOTO (resume from last photo)" >> "$LOG_FILE"
    START=$(date "+%Y-%m-%d %H:%M:%S")
    ./re "$PHOTO" > re.out 2>&1
    RE=$?
    END=$(date "+%Y-%m-%d %H:%M:%S")
    echo "  re start=$START end=$END exit=$RE (see $dir/re.out)" >> "$LOG_FILE"
    cd "$BASE_DIR"
done

echo "" >> "$LOG_FILE"
echo "=== photo verification ===" >> "$LOG_FILE"
FAIL=0
for dir in "${DIRS[@]}"; do
    [ -d "$BASE_DIR/$dir" ] || continue
    np=$(photo_count "$BASE_DIR/$dir")
    if [ -f "$BASE_DIR/$dir/restart_photo" ] && [ -s "$BASE_DIR/$dir/restart_photo" ]; then
        rp=yes
    else
        rp=no
    fi
    if [ "$np" -eq 0 ]; then
        FAIL=$((FAIL+1))
        echo "[FAIL] $dir nonzero_photos=$np restart_photo=$rp" >> "$LOG_FILE"
    else
        echo "[OK]   $dir nonzero_photos=$np restart_photo=$rp" >> "$LOG_FILE"
    fi
done
echo "" >> "$LOG_FILE"
echo "photos check failures: $FAIL" >> "$LOG_FILE"
echo "Done at $(date)" >> "$LOG_FILE"
