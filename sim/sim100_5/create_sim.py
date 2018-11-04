#!/usr/bin/env python
#
#
import os

for salpha in range(4,11):
    alpha = salpha / 10.
    os.system("cp parmfile.simtemplate parmfilesim_{0}".format(alpha))
    os.system("perl -p -i -e 's/MLFALPHA/{0}/' parmfilesim_{0}".format(alpha))
    os.system("cat parmfilesim_{0}  | ~/src/simulator/migtree2 | ~/src/simulator/migdata; cp infile infile.{0}".format(alpha))
    os.system("cp parmfile.template parmfile_{0}".format(alpha))
    os.system("perl -p -i -e 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))

    
