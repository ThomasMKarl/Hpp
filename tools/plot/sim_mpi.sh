#!/bin/bash
echo "backup..."
rm -f *~
mv -f ising1.dat ising1.dat.bak
mv -f ising4.dat ising4.dat.bak
mv -f ising9.dat ising9.dat.bak 
echo "done."

echo "simulating..."
for i in 6 12 18 30 42 60 90 138 204; do
    echo " "
    echo "grid size $i"
    mpirun --hostfile hostfile -np 1 ./ising $i $i hot_start  10 1 0 10000 10 >> ising1.dat
    mpirun --hostfile hostfile -np 4 ./ising $i $i hot_start  10 1 0 10000 10 >> ising4.dat
    mpirun --hostfile hostfile -np 9 ./ising $i $i hot_start  10 1 0 10000 10 >> ising9.dat
done
echo "done."
