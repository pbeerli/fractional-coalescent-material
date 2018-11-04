#!/usr/bin/env python3
# simulator for genealogies of different nestsite quality
#
# (c) peter beerli
#
# uses python 3 syntax but should run on 2
#
#
#
# usage: nestsite.py [-h] [-v] [-m MUTATIONRATE] [-N POPSIZE] [-n SAMPLESIZE]
#                    [-l LOCI] [-s SITES] [-g GENERATIONS] [-a ALPHA]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -v, --verbosity       increase output verbosity
#   -m MUTATIONRATE, --mutationrate MUTATIONRATE
#                         specify the mutation rate to scale the tree
#   -N POPSIZE, --populationsize POPSIZE
#                         specify the number of individuals in the population
#   -n SAMPLESIZE, --samplesize SAMPLESIZE
#                         specify the number of individuals in the sample
#   -l LOCI, --loci LOCI  specify how many sites should be asked for to simulate
#   -s SITES, --sites SITES
#                         specify how many sites should be asked for to simulate
#   -g GENERATIONS, --generations GENERATIONS
#                         number of generations to run
#   -a ALPHA, --heterogeneity ALPHA
#                         alpha indicates the relative success of enviroment for
#                         breeding, for example [1 2 3 10]
#
# This will print a tree with commented instructions for a companion program migdata like this:
# #SN
# #1365906878
# #1
# #10
# #1 1000 2.000000
# # rate among sites for locus 0 (1.000000)
# #=
# #$ locus 1
# #$ 0.0001
# (((a300005:0.002,(a300019:0.0013,a300023:0.0013):0.0007):0.0003,a300038:0.0023):0.0085,(a300047:0.0034,(a300075:0.0011,(a300083:0.0005,(a300088:0.0003,a300092:0.0003,a300099:0.0003):0.0002):0.0006):0.0023):0.0074):0.0;
# 
# an example: python nestsite.py > tree ; cat tree | migdata
# will then generat a sequence data file with 10 individuals using the genealogy in Newick format,
# with mutation rate 0.0001
# 
from __future__ import print_function
from __future__ import division

import numpy as np
import scipy as sc
import sys
import os
import time

import tree as Tree

# generate first generation
def startup(n):
    tmp =  [Tree.Node(i) for i in range(n)]
    for i in range(n):
        tmp[i].age = 0
    return tmp

# alpha is a vector: it has the size of number of site types
# n is the number of individuals in the population
# returns a weight for each individual based on dirichlet alpha
def nestsitequality(alpha,n):
    weight =  np.random.dirichlet(alpha)
    #eprint("Weights:", weight)
    newweight =  np.random.choice(weight,size=n,replace=True)
    sumweight = np.sum(newweight)
    return newweight/sumweight

# there are n individuals each has a list of ancestors
# each individual draws from the parent generation, this will be a list
# of parentIDs these will be then matched to the last individual on
# record (the parent), and if there are multiple they will be
# linked using a sub array
def next_generation(individuals, npop,weight):
    x = [i.name for i in individuals]
    return np.random.choice(x, npop, p=weight,replace=True)

def hookup(f2names, f1):
    global counter
    newf1 = []
    for fi in f2names:
        ii = 0
        while fi != f1[ii].name: 
            ii += 1
        x = Tree.Node(counter)
        counter += 1
        x.set_parent(f1[ii])
        f1[ii].set_daughter(x)
        #x.add_age(1)
        #print(ii,f1[ii])
        newf1.append(x)
    return newf1

# simulates the genealogy of a whole population then selects n individuals
# at random
#
def sim(popsize, n, alphas, generations, mutationrate):
    alpha = np.array(alphas)
    # out = open('testtree','w+')
    weight = nestsitequality(alpha,popsize)
    #print("nest site quality:",weight)
    f1 = startup(popsize)
    starttime=time.time()
    for g in range(generations):
        #eprint("gen: {}".format(g))
        weight = nestsitequality(alpha,popsize)
        f2names = next_generation(f1,popsize,weight)
        f2names.sort()
        f1 = hookup(f2names,f1)
        for fi in f1:
            fi.age = g
            #print("#%",fi.name,fi.age,fi.ancestor.name)
    endtime = time.time()
    oldroot = Tree.lastancestor(f1[0])
    for i in range(1,popsize+1):
        #eprint("checked {} of {}".format(i,popsize))
        root = Tree.lastancestor(f1[i-1])
        if root != oldroot:
            eprint("more than one root with contemporary tips");
            print("more than one root with contemporary tips");
    eprint("#runtime: found {} roots".format((root!=oldroot) + 1) ,endtime-starttime)
    starttime = time.time()
    tree = Tree.Tree(root)
    # this helps to prune the large popuulation tree
    countn=0
    for fi in np.random.choice(f1, n,replace=False):
        fi.deadend = 0
        countn+=1
        #eprint("#",fi.age,fi.name, end=" ")
    eprint("N=",len(f1)," n=",countn)
    tree.fixdeadends_norec(f1,root)
    endtime = time.time()
    eprint("#runtime: fixdeadends",endtime-starttime)        
       
    #tree.prune(root)
    starttime = time.time()
    tree.prune_norec(f1, root)
    endtime = time.time()
    eprint("#runtime",endtime-starttime)        

    #for fi in f1:
    #    print("#2 ",fi.age,fi.name)
    tree.set_blength(root)
    #tree.set_blength_norec(f1,root)
    tree.myprintprune(root,mutationrate)
    
def eprint(*args, **kwargs):
    global verbose
    if verbose:
        print(*args, file=sys.stderr, **kwargs)        

def eprint2(*args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)        
        sys.stderr.flush()

def print_header(n,loci,sites):
    print("#SN")
    print("#{}".format(np.random.randint(1,2**31)))
    print("#1")
    print("#{}".format(n)) 
    print("#{} {} 2.000000".format(loci,sites))
    for locus in range(loci):
        print("# rate among sites for locus 0 (1.000000)")
        print("#=")

def print_treeheader(locus, mutation):
    print("#$ locus {}".format(locus))
    print("#$ {}".format(mutation))

    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbosity", action="count",default=0,
                    help="increase output verbosity")
    parser.add_argument('-m','--mutationrate', dest='mutationrate', type=float,default=0.0001,
                        help='specify the mutation rate to scale the tree')
    parser.add_argument('-N','--populationsize', dest='popsize', type=int,default=100,
                        help='specify the number of individuals in the population')
    parser.add_argument('-n','--samplesize', dest='samplesize', type=int,default=10,
                        help='specify the number of individuals in the sample')
    parser.add_argument('-l','--loci', dest='loci', type=int,default=1,
                        help='specify how many sites should be asked for to simulate')
    parser.add_argument('-s', '--sites', dest='sites', type=int,default=1000,
                        help='specify how many sites should be asked for to simulate')
    parser.add_argument('-g', '--generations', dest='generations', type=int,default=3000,
                        help='number of generations to run')

    parser.add_argument('-a', '--heterogeneity', dest='alpha',default='1 2 4',
                        help='alpha indicates the relative success of enviroment for breeding, for example [1 2 3 10]')

    args = parser.parse_args()
    verbose = args.verbosity
    Tree.verbose = verbose
    popsize = args.popsize
    mutationrate = args.mutationrate
    n = args.samplesize
    loci = args.loci
    sites = args.sites
    generations = args.generations
    alphas = list(map(float,args.alpha.split()))
    
    eprint("OPTIONS:")
    eprint("popsize=",popsize)
    eprint("mutationrate=",mutationrate)
    eprint("n=",n)
    eprint("verbose=",'yes')
    eprint("loci=",loci)
    eprint("sites=",sites)
    eprint("generations=",generations)
    eprint("alphas=",alphas)
    
    #sys.exit()
    sys.setrecursionlimit(2500)
    counter = popsize
    print_header(n,loci,sites)
    for locus in range(1,loci+1):
        print_treeheader(locus,mutationrate)
        sim(popsize, n, alphas, generations, mutationrate)
        print(":0.0;")       
        eprint2(locus,end=" ")
    eprint2("")    
    eprint("OPTIONS:")
    eprint("popsize=",popsize)
    eprint("mutationrate=",mutationrate)
    eprint("n=",n)
    eprint("verbose=",'yes')
    eprint("loci=",loci)
    eprint("sites=",sites)
    eprint("generations=",generations)
    eprint("alphas=",alphas)
