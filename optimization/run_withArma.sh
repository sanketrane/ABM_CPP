#!/bin/sh

while getopts m: flag
do
    case "${flag}" in
        m) filename=${OPTARG};;
    esac
done

c++ -std=c++17 -larmadillo ${filename}.cpp -o ${filename} 
./${filename} > ${filename}_res.out


echo "---END CODE---"
