#!/bin/sh
./exptst << eof 
100
$1
$2
5
n
n
1 2
eof
mv pgplot.ps exptst.ps
gv exptst.ps
