#!/usr/bin/env python

## altuntas@ucar.edu 2019
descr = """
extract xmlchange commands from a caseroot
"""

import os
import argparse
import math
from collections import namedtuple
from operator import attrgetter, itemgetter


###############################################################################

Routine = namedtuple("Routine", "index total self perc child name line")

def read_amplxe_results(results_path):

    with open(results_path,'r') as results:

        # read header:
        results.readline()

        # read column lengths:
        line = results.readline()
        cw = [len(l)+2 for l in line.split()]
        slices = []
        for i in range(5):
            slices.append( slice( sum(cw[:i]), sum(cw[:i+1]) ) )

        routines = dict()

        # now read the entries
        i=2
        exetime = 0.0
        for line in results:
            i += 1 # line number


            if len(line.split())==0:
                continue # empty line
            elif line.split()[0]=="Index":
                break # end of stats section

            # index
            index = str(line[slices[0]])

            try:
                # total time
                try:
                    total = float(line[slices[1]])
                except:
                    total = 0

                # self time
                self = float(line[slices[2]])
            
                # child time
                child = float(line[slices[3]])

                # name
                name = line[slices[4]].strip()

                if index.strip()=="[1]":
                    exetime = (self+child)/total*100.0
            except:
                continue # irregular line

            ## save only indexed entries
            if index.strip() == "":
                continue

            routines[name] = Routine(index=index, total="{:7.3f}".format(total),
                                                  self="{:10.4f}".format(self),
                                                  perc="{:10.4f}".format(self/exetime*100),
                                                  child="{:10.4f}".format(child),
                                                  name=name, line=i)
    return routines

def print_most_expensive(routines, n, args):

    # from dict to list:
    routines_list = [routines[r] for r in routines]
    routines_list_sorted = sorted(routines_list, key=attrgetter('self'), reverse=True)

    print("{:8} {:10s} {:10s} {:10s} {:11s} {:10s}".\
            format("index"," total%", "  self%", " self", " child", "name"))

    nprinted = 0
    for i in range(len(routines_list_sorted)):
        r = routines_list_sorted[i]
        if args.t>float(r.perc):
            continue
        print("{:8} {:10s} {:10s} {:10s} {:11s} {}".\
            format(r.index,r.total,r.perc, r.self, r.child, r.name))
        nprinted +=1
        if nprinted==args.n:
            break
        

def compare_results(routines1, routines2, n, args):

    routines_diff = []
    routines = set(routines1.keys()).union(set(routines2.keys()))
    for rname in routines:
        if rname in routines1:
            r1_self = float(routines1[rname].self)
            r1_child = float(routines1[rname].child)
            r1_selfperc = float(routines1[rname].perc)
        else:
            r1_self = math.nan
            r1_child = math.nan
            r1_selfperc = math.nan
        if rname in routines2:
            r2_self = float(routines2[rname].self)
            r2_child = float(routines2[rname].child)
            r2_selfperc = float(routines2[rname].perc)
        else:
            r2_self = math.nan
            r2_child = math.nan
            r2_selfperc = math.nan

        if not (math.isnan(r1_self) or math.isnan(r2_self)):
            if r2_self != 0:
                ratio = r1_self/r2_self
            else:
                ratio = math.nan
        else:
            ratio = math.nan

        routines_diff.append([
            rname, 
            ratio,
            r1_self,
            r2_self,
            r1_child,
            r2_child,
            (r1_self+r1_child)/(r2_self+r2_child),
            r1_selfperc,
            r2_selfperc])

    if args.s == "self":
        ratio_ix = 1
    elif  args.s == "total":
        ratio_ix = 6
    else:
        raise RuntimeError("Invalid -s option: specify 'self' or 'total' ")
    routines_diff_sorted = sorted([x for x in routines_diff if not math.isnan(x[ratio_ix])], key=itemgetter(ratio_ix), reverse=True) + [x for x in routines_diff if math.isnan(x[ratio_ix])]
    print("{:11s} {:11s} {:9} {:9} {:9} {:10s} {}".\
            format(" ratio", " self1", " self2", "child1", "child2", "self1%", "self2%", "name"))
    nprinted = 0
    for i in range(len(routines_diff_sorted)):
        rname = routines_diff_sorted[i][0]
        ratio = routines_diff_sorted[i][ratio_ix]
        self1 = routines_diff_sorted[i][2]
        self2 = routines_diff_sorted[i][3]
        child1 = routines_diff_sorted[i][4]
        child2 = routines_diff_sorted[i][5]
        self1perc = routines_diff_sorted[i][7]
        self2perc = routines_diff_sorted[i][8]

        if not ( math.isnan(self1perc) or math.isnan(self2perc) ):
            if args.t>self1perc and args.t>self2perc:
                continue
        elif not math.isnan(self1perc):
            if args.t/2>self1perc:
                continue
        elif not math.isnan(self2perc):
            if args.t/2>self2perc:
                continue

        print("{:10f} {:9} {:9} {:9} {:9} {:9} {:9} {}".\
            format(ratio,self1,self2,child1,child2,self1perc,self2perc,rname))
        nprinted +=1
        if nprinted==args.n:
            break

def print_all_indices(routines, args):
    # from dict to list:
    routines_list = [routines[r] for r in routines]
    routines_list_sorted = sorted(routines_list, key=attrgetter('index'))

    print("{:8} {:10s} {:10s} {:10s} {:11s} {:10s}".\
            format("index"," total%", "  self%", " self", " child", "name"))

    for i in range(len(routines_list_sorted)):
        r = routines_list_sorted[i]
        if args.t>float(r.perc):
            continue
        print("{:8} {:10s} {:10s} {:10s} {:11s} {}".\
            format(r.index,r.total,r.perc, r.self, r.child, r.name))


def main():

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('-f1', metavar='', type=str, required=True,
                        help='results_txt file of a run')
    parser.add_argument('-f2', metavar='', type=str, required=False,
                        help='results_txt file to compare with f1')
    parser.add_argument('-t', metavar='', type=float, required=False,
                        help='threshold for self percentage of f1', default=0.0)
    parser.add_argument('-n', metavar='', type=int, required=False, default=20,
                        help='number of routines to print')
    parser.add_argument('-s', metavar='', type=str, required=False, default="self",
                        help='sort by: "self" or "total" ')
    args = parser.parse_args()

    routines1 = read_amplxe_results(args.f1)
    if args.f2:
        # two results file provided
        routines2 = read_amplxe_results(args.f2)
        compare_results(routines1, routines2, args.n, args)
    else:
        # single results file provided
        print_most_expensive(routines1, args.n, args)

if __name__ == "__main__":
    main()

