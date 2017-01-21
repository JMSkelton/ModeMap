#!/bin/bash

# Wrapping the loop body in a function saves typing (!).

function Run()
{
  temp=${1}

  # This script can take a long time to run -> print a "status message" to show that something is happening.

  echo "Running: T = ${temp} K"

  # Modify input.dat to overwrite the temperature.
  # Note that this code assumes the layout of the example input file (i.e. with the temperature appearing on line 22).

  head -21 "input.dat" > "input_temp.dat"
  echo "${temp}" >> "input_temp.dat"
  tail +23 "input.dat" >> "input_temp.dat"

  mv "input_temp.dat" "input.dat"

  # Run the solver.

  ./schrodinger-solve1d-anharmonic.exe > out.dat

  # Extract the effective (renormalised) harmonic frequency.

  freq=`grep "Omega" "out.dat" | tail -3 | head -1 | awk '{print $3}'`

  # Append the temperature and effective frequency to "t-vs-omega.dat".

  echo "${temp}  ${freq}" >> "t-vs-omega.dat"

  # Merge the output into "full-out.dat"

  echo "Temperature: ${temp}" >> "full-out.dat"
  echo "" > "full-out.dat"

  cat "out.dat" >> "full-out.dat"

  echo "" >> "full-out.dat"
}

# Set up an output file to hold the calculated effective frequencies.

echo "# T [K]  Frequency [THz]" > "t-vs-omega.dat"

# We want to calculate the frequencies at low temperature over a fine temperature range...

for (( i=0; i<10; i++ ))
do
  Run ${i}
done

# ... but the frequencies at higher temperatures can be calculated over a coarser range.

for (( i=1; i<=100; i++ ))
do
  Run $(($i * 10))
done
