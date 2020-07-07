#!/bin/bash
echo "backup..."
rm -f *~
mv -f ising.dat ising.dat.bak
echo "done."

echo "simulating..."
echo "#Temperature    Energy    Error" > ising2.dat
for t in 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 6.25 6.5 6.75 7.0 7.25 7.5 7.75 8.0; do
    echo " "
    echo "Temperature: $t x J/k_b"
    ./ising $t 1.0 0.0 30 30 30 1000 "cold_start" >> ising2.dat
done
echo "done."
