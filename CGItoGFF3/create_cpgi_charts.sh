#!/bin/bash

REPORTS_DIR=/home/brock/bionformatics/arabidopsis/TAIR10/reports
IMAGES_DIR=./
REPORTS_EXT=cpgreport

cat $REPORTS_DIR/*.$REPORTS_EXT | ./CGItoGFF3 > $REPORTS_DIR/islands.gff3

for name in Chr1 Chr2 Chr3 Chr4 Chr5 mitochondria chloroplast
do
    gt sketch -force -width 2400 -seqid $name $name-islands.png $REPORTS_DIR/islands.gff3
done
