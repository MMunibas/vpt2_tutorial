#!/bin/bash
for i in *.inp
do
    filename="${i%.*}"
    if [ ! -f $filename.out ]
    then
        echo "molprosub" $filename.inp
    fi
done
