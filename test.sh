#!/bin/bash
set -e

NPLANETS=${1}
TIMESTEPS=${2}

CORRECT_ICON="\U2705"
WRONG_ICON="\U274C"

./main.exe "$NPLANETS" "$TIMESTEPS" > "output/$NPLANETS-$TIMESTEPS.txt"

# Extract the final location from the output file
final_location=$(tail -n 1 "output/$NPLANETS-$TIMESTEPS.txt" | awk '{print $10, $11}')

# Extract the final location from the correct output file
correct_final_location=$(tail -n 1 "output/${NPLANETS}-${TIMESTEPS}_correct.txt" | awk '{print $10, $11}')

if [ "$final_location" == "$correct_final_location" ]; then
    echo -e "$CORRECT_ICON $NPLANETS-$TIMESTEPS passed"
else
    echo -e "$WRONG_ICON $NPLANETS-$TIMESTEPS failed"
fi