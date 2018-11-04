#!/usr/bin/env python
#
#
import os

for salpha in range(4,10):
    alpha = salpha / 10. + 0.05
    os.system("cp parmfile.realtemplate parmfile_{0}".format(alpha))
    os.system("perl -p -i -e 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))
    #os.system("~/src/migrate-mlf/migrate-n  parmfile_{0} -nomenu".format(alpha))
    os.system("submit ./migrate-mlf-mpi parmfile_{0} log_{0} 41 4:00:00".format(alpha))
for salpha in range(4,11):
    alpha = salpha / 10.
    os.system("cp parmfile.realtemplate parmfile_{0}".format(alpha))
    os.system("perl -p -i -e 's/MLFALPHA/{0}/' parmfile_{0}".format(alpha))
    #os.system("~/src/migrate-mlf/migrate-n  parmfile_{0} -nomenu".format(alpha))
    os.system("submit ./migrate-mlf-mpi parmfile_{0} log_{0} 41 4:00:00".format(alpha))

    
