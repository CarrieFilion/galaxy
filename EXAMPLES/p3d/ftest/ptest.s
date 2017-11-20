#!/bin/sh
rm -f run$1.lis
ptest << eof
$1
eof
cat run$1.lis
