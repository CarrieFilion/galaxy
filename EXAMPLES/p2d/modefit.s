#!/bin/sh
# finds the best fit mode to the data output from run100
modefit << eof
100
3
n
n
1 2
lgsp
2
n
n
21 56
10
161
0
n
1
.1
n
n
n
n
2
2
.1
n
n
y
y
n
n
0
n
eof
gv pgplot.ps
