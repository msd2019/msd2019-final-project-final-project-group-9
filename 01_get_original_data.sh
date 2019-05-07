#!/bin/bash
FP="data/original"

#Downloads the data into a zip file
Curl -L0 http://tuvalu.santafe.edu/~aaronc/facultyhiring/replicationData_all.zip > "$FP/replicationData_all.zip"

#unzip the data
unzip "$FP/replicationData_all.zip" -d $FP


