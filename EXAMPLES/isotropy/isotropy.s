#!/bin/sh
isotropy << eof 
$1
100
4
y
18 12
n
2 1
eof
mv pgplot.ps run$1.ps
gv run$1.ps
