#!/bin/bash
for i in *.out
do
    normal_termination=$(grep -c "Molpro calculation terminated" $i) 
    if [ $normal_termination -eq 0 ]; then
       echo $i
       #rm $i
    fi
done
