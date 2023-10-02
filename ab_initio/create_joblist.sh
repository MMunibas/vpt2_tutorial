#!/bin/bash
for i in *.inp
do
    filename="${i%.*}"
    if [ ! -f $filename.out ]
    then
        echo "molpro2022_par" $filename.inp
    fi
done
