#/usr/bin/env python
#
'''
SLURM functionality 
'''
import os

def is_inqueue(file):
    '''
    Checks whether a particular SLUMR batchfile is running, file should be a path:
    for example: "x0dx10/log_0_3.25_10.slurm"
    '''
    a = int(os.popen('squeue -u beerli -o "%all" | grep '+file+' | wc -l').read())
    if a > 0:
        return True
    else:
        return False

if __name__ == '__main__':
    import sys
    print is_inqueue(sys.argv[1])

