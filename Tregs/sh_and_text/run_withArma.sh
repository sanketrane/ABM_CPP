#!/bin/sh

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

c++ -std=c++17 Models/${modelname}.cpp -o ${modelname} -larmadillo
./${modelname} m1 > ${modelname}_res.out


echo "---END CODE---"
