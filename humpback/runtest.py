#!/usr/bin/env python
#
#
import os

for salpha in range(4,10):
    alpha = salpha / 10.
    os.system("cp parmfile.realtemplate parmfile_{0}".format(alpha))
    os.system("sed -i'' 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))
    os.system("submit ./migrate-mlf-mpi parmfile_{0} log_{0} 11 4:00:00 backfill2".format(alpha))
    alpha = alpha + 0.05
    os.system("cp parmfile.realtemplate parmfile_{0}".format(alpha))
    os.system("sed -i'' 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))
    os.system("submit ./migrate-mlf-mpi parmfile_{0} log_{0} 11 4:00:00 backfill2".format(alpha))

alpha = 1.0
os.system("cp parmfile.realtemplate parmfile_{0}".format(alpha))
os.system("sed -i'' 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))
os.system("submit ./migrate-mlf-mpi parmfile_{0} log_{0} 11 4:00:00 backfill2".format(alpha))

