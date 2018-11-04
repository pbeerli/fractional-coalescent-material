#!/usr/bin/env python
#
#
import os
import slurm
#for salpha in [5,6,7,8,9,10]:
os.system("mpirun -np 9 ~/src/migrate-mlf/migrate-n-mpi parmfile_1.0_1.0s -nomenu")
#os.system("launch ./migrate-mlf-mpi  parmfile_1.0_1.0s log_1.0_1.0s 21 40:00:00"
for salpha in [10]:
    alpha = salpha / 10.
    for salpha2 in [10,9]:
        alpha2 = salpha2 / 10.
        outfile = "outfile_{0}_{1}.pdf".format(alpha,alpha2)
        log = "log_{0}_{1}".format(alpha,alpha2)
        mydir = os.getcwd()
        if not os.path.isfile(outfile):
            if not slurm.is_inqueue(mydir+"/"+log):
                os.system("cp parmfile.template parmfile_{0}_{1}".format(alpha,alpha2))
                os.system("perl -p -i -e 's/MLFALPHA/{0}/' parmfile_{0}_{1}".format(alpha,alpha2))
                os.system("perl -p -i -e  's/MLF2ALPHA/{1}/' parmfile_{0}_{1}".format(alpha,alpha2))
                #os.system("submit ./migrate-mlf-mpi  parmfile_{0}_{1} log_{0}_{1} 21 40:00:00".format(alpha,alpha2))
                os.system("mpirun -np 4 ~/src/migrate-mlf/migrate-n-mpi parmfile_{0}_{1} -nomenu".format(alpha,alpha2))
        

    
