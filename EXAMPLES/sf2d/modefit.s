#!/bin/sh
modefit << eof
100    # run number
3
n
y
1 2
lgsp   # select log spirals
2      # sectoral harmonic
n      # do not force zero growth rate
n      # do not force zero pattern speed
21 56  # range of log spirals
5      # first moment
55     # analysis steps
2      # draw selected data
n      # do not reselect
1      # number of modes to be fitted
.1     # exponential scale factor
y      # change initial guess
.6 .2  # initial guess
n      # show this fit
n      # show mode shapes
n      # do not save this mode shape
n      # do not save this fit
2      # change number of modes
2      # 2 modes
.1     # exponential scale factor
n      # do not change initial guess
n      # do not change initial guess
y      # show this fit
y      # show mode shapes
n      # do not save this mode shape
n      # do not save this fit
0      # done
n      # no other data
eof
mv pgplot.ps modes.ps
gv modes.ps
