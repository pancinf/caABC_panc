#!/bin/bash
#This script runs all download scripts

./downloadRaw.sh > ../../data/raw/downloadRaw.log 2>&1
./downloadRef.sh > ../../data/ref/downloadRef.log 2>&1
./downloadGTEx.sh > ../../data/GTEx/raw/downloadGTEx.log 2>&1
./downloadCage.sh > ../../data/cage/signals/downloadCage.log 2>&1