#!/bin/sh

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

c++ -std=c++17 ${modelname}.cpp -o ${modelname}
./${modelname} > ${modelname}_res.out

echo "---END CODE---"