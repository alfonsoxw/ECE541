from __future__ import print_function
import random
import pdb
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(u'-seed', metavar  = u'initial seed', dest=u'seed', required=True)
parser.add_argument(u'-steps', metavar = u'steps', dest=u'steps', required=True)
parser.add_argument(u'-die', metavar = u'die', dest=u'die', required=True)
parser.add_argument(u'-population', metavar = u'population', dest=u'pop', required=True)

args = parser.parse_args()

random.seed(args.seed)

Steps = int( args.steps )

pDie          = float(args.die)
pNoChildren   = (1.0-pDie)/4
pOneChild     = (1.0-pDie)/4
pTwoChildren  = (1.0-pDie)/4
pThree        = (1.0-pDie)/4

cdf    = [0]*5
cdf[0] = pDie
cdf[1] = cdf[0]+pNoChildren
cdf[2] = cdf[1]+pOneChild
cdf[3] = cdf[2]+pTwoChildren
cdf[4] = 1.0

mean =  pDie*0
mean += pNoChildren*1
mean += pOneChild*2
mean += pTwoChildren*3
mean += pThree*4

print('mean individual update', mean)

def sample():
    u = random.random()
    idx = 0 
    while cdf[idx] < u:
        idx += 1

    return idx

X = int(args.pop)

for s in xrange(0,Steps):
    if X <= 0:
        break
    print('step',s,'population',X)
    Y=0
    for ydx in xrange(0, X):
        Y += sample()
    X = Y 
    

