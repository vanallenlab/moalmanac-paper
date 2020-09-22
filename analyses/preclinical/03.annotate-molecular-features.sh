#!/bin/bash

MOALMANAC_DIR=$1
MOALMANAC_VENV=$2
WD=$PWD

cp almanac-gdsc-mappings.json $1/datasources/preclinical/almanac-gdsc-mappings.json
cp formatted/cell-lines.summary.txt $1/datasources/preclinical/cell-lines.summary.txt
cp formatted/sanger.gdsc.txt $1/datasources/preclinical/sanger.gdsc.txt

cp annotate-variants.py annotate-copy-numbers.py annotate-fusions.py $MOALMANAC_DIR
cd $MOALMANAC_DIR

$MOALMANAC_VENV annotate-variants.py --directory $WD
$MOALMANAC_VENV annotate-copy-numbers.py --directory $WD
$MOALMANAC_VENV annotate-fusions.py --directory $WD

rm annotate-variants.py annotate-copy-numbers.py annotate-fusions.py
