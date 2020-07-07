#!/bin/bash
echo "backup..."
rm -f *~
mv -f time.dat time.dat.bak
echo "done."

echo "simulating..."
echo "#Temperature    Energy    Error" > time.dat
for g in 12 24 36 48 60 72 84 96 108 120; do
    echo " "
    echo "Size: $g^3"
    ./ising 2.0 1.0 0.0 $g $g $g 1000 "cold_start" >> time.dat
done
echo "done."
