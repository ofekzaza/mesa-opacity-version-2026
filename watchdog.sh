#!/usr/bin/env bash
# watchdog.sh [run_dir] [idle_minutes]
# Runs ./rn (first time) or ./re (restart). If ./star stops writing to LOGS/
# for idle_minutes while still running, kill it and run ./re again.
set -u
cd "${1:-.}" || exit 1
idle=$(( ${2:-60} * 60 ))
trap 'pkill -f "./star"; exit 130' INT TERM

while :; do
  if [ -n "$(ls photos 2>/dev/null | head -1)" ]; then
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
