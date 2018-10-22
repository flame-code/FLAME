#!/bin/bash

FOLDERS=$(ls POSCAR-*)

for i in $FOLDERS 
do
echo $i
mkdir DIR_$i
cd DIR_$i
cp ../SAMPLE/* .
cp ../$i .
#cp $i POSCAR
sed -i '5 a Ca F' POSCAR*
poscar2yaml.py $i posinp.yaml
#echo "POSCAR"> in
#echo "T">> in
#echo "Na Cl" >> in
#POSCAR2ascii.x < in
#poscar2ascii.py POSCAR > POSCAR.ascii
#ascii2acf.py POSCAR.ascii  > posinp.acf
../../../../flame >&o1
cd ..
echo '---------------------------------------------------------------------'
done


