#!/usr/bin/env python

from plot_utils import popSizeStepPlot

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-I",       "--msmcFinal",            help="Specify MSMC result file",                      required = True,            metavar="MSMC result file")
parser.add_argument("-u",       "--mutationRate",         help="Specify mutation rate",                         required = True,            metavar="Mutation rate")
parser.add_argument("-g",       "--generations",          help="Specify generations per year",                  required = True,            metavar="generations")
args = parser.parse_args()

msmcFile        = args.msmcFinal
mut             = float(args.mutationRate)
gener           = float(args.generations)

# x contains the left point of each step-segment in years
# y contains the effective population size

(x,y) = popSizeStepPlot(msmcFile, mu=mut, gen=gener)

print('leftX\tpopSize')
for i in range(len(x)):
    print(x[i], y[i])
