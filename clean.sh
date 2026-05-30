#!/bin/bash

for d in */; do
    # Check if it is a directory
    if [ -d "$d" ]; then
        # Check if clean exists and is executable
        if [ -x "${d}clean" ]; then
            echo "Running ./clean in ${d}"
            (cd "$d" && ./clean)
        # Fallback if clean exists but is not executable
        elif [ -f "${d}clean" ]; then
            echo "Running bash clean in ${d}"
            (cd "$d" && bash clean)
        fi
    fi
done
