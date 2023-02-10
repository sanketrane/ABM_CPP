#!/bin/sh

c++ -std=c++17 mod_trial.cpp -o mod_trial
./mod_trial > mod_res.out

echo "---END CODE---"