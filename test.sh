#!/bin/bash
set -e

NPLANETS=${1}
TIMESTEPS=${2}

CORRECT_ICON="\U2705"
WRONG_ICON="\U274C"

# Run sim and save output
./main.exe "$NPLANETS" "$TIMESTEPS" > "output/$NPLANETS-$TIMESTEPS.txt"

# Compare with correct output
if ./compare_outputs "output/$NPLANETS-$TIMESTEPS.txt" "output/${NPLANETS}-${TIMESTEPS}_correct.txt"; then
    echo -e "$CORRECT_ICON $NPLANETS-$TIMESTEPS passed"
else
    echo -e "$WRONG_ICON $NPLANETS-$TIMESTEPS failed"
fi