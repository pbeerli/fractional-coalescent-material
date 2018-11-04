import sys
import os
alpha = sys.argv[1]
fo = sys.argv[2]
#os.system("rm gugus")
for j in [1,2,3,4,5]:
    for i in range(1,10):
        f = "grep ^'      {}     '  {}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(i, fo, alpha)
        print "{:4d} {}".format(i,f)
        os.system(f)
    for i in range(10,100):
        f = "grep ^'     {}     '  {}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(i, fo, alpha)
        print "{:4d} {}".format(i,f)
        os.system(f)

    f = "grep ^'    {}     '  {}/outfile_{}_*  | bf.py  | grep  '   0.00  ' >>  gugus".format(100, fo, alpha)
    print "{:4d} {}".format(100,f)
    os.system(f)

