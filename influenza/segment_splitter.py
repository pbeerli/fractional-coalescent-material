#!/usr/bin/env python
#
import os
import sys
for segment in range(1,9):
    os.system("grep -A1 Segment:{0} influenza_h1n1_mexico_2014.sorted | grep -v '\v\v' > tmp ;  mafft --auto tmp > tmp2; fasta2mig tmp2 > influenza_h1n1_mexico_2014_segment{0}".format(segment))
    
    
    
