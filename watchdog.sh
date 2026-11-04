#!/usr/bin/env bash
# watchdog.sh [--re] [run_dir] [idle_minutes]
# Runs ./rn (first time) or ./re (restart). If ./star stops writing to LOGS/
# for idle_minutes while still running, kill it and run ./re again.
#   --re   always start with ./re, never ./rn (resume from the last photo)
set -u
force_re=0
args=()
for a in "$@"; do
  case "$a" in
    --re|-re) force_re=1 ;;
    *) args+=("$a") ;;
  esac
done
cd "${args[0]:-.}" || exit 1
idle=$(( ${args[1]:-60} * 60 ))
trap 'pkill -f "./star"; exit 130' INT TERM

while :; do
  if [ "$force_re" = 1 ] || [ -n "$(ls photos 2>/dev/null | head -1)" ]; then
    ./re >>watchdog.log 2>&1 &
  else
    ./rn >>watchdog.log 2>&1 &
  fi
  w=$!
  killed=0
  while kill -0 "$w" 2>/dev/null; do
    sleep 30
    kill -0 "$w" 2>/dev/null || break
    f=$(ls -t LOGS 2>/dev/null | head -1)
    [ -n "$f" ] || continue
    m=$(stat -f %m "LOGS/$f" 2>/dev/null) || continue
    if [ $(( $(date +%s) - m )) -ge "$idle" ]; then
      echo "watchdog: no output for $idle s -> killing star, restarting from photo"
      pkill -f './star'
      killed=1
      break
    fi
  done
  wait "$w" 2>/dev/null
  [ "$killed" = 1 ] || break
done
