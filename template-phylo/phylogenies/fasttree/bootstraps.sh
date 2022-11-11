#!/bin/bash
DIR=$(dirname $2)
FA=$1
OUT=$2
cd $DIR
cp $FA bs.fa
raxmlHPC-AVX -f j -b 38491 -# 100 -s bs.fa -n bs -m GTRCAT
cat bs.fa.BS* >$OUT
