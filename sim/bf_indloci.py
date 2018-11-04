import sys
import os
alpha = sys.argv[1]
#os.system("rm gugus")
for j in [1,2,3,4,5]:
    for i in range(1,10):
        f = "grep ^'      {}     '  sim100_{}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(i, j, alpha)
        print "{:4d} {}".format(i,f)
        os.system(f)
    for i in range(10,100):
        f = "grep ^'     {}     '  sim100_{}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(i, j, alpha)
        print "{:4d} {}".format(i,f)
        os.system(f)

    f = "grep ^'    {}     '  sim100_{}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(100, j, alpha)
    print "{:4d} {}".format(100,f)
    os.system(f)

