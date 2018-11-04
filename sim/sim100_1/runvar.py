#!/usr/bin/env python
import os
for i in range(1,11):
    a = i/10.
    os.system("python /Users/beerli/src/migrate-new/countvarsites.py infile."+str(a)+" | grep All | awk '{print "+str(a)+",$0}' >> variable_sites")
