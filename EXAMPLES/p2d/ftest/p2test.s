#!/bin/sh
ptest << eof 
$1
eof
cat run$1.lis
