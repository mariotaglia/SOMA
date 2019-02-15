#!/bin/bash

FLAGS="-linux --brace-indent 4 -bli4 -npcs -ppi0 --line-length 120 -sob --indent-level 4 -bl -nce -nut --cuddle-do-while"

#HEADER files
for file in *.h
do
    indent $FLAGS $file
done

#SOURCE files
for file in *.c
do
    indent $FLAGS $file
done

#EXTRA files
for file in soma_config.h.in
do
    indent $FLAGS $file
done
    
