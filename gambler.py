#from __future__ import print_function
import pdb
import math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(u'-pwin', metavar  = u'probability of winning', dest=u'pwin', required=True)
parser.add_argument(u'-jackpot', metavar = u'jackpot', dest=u'jackpot', required=True)
parser.add_argument(u'-stake', metavar = u'stake', dest=u'stake', required=True)

args = parser.parse_args()

pUp = float(args.pwin)
pDn = 1-pUp

Jackpot  = int(args.jackpot)
Stake = int(args.stake) 

def nxtStep(pi):
    nxtPi = [0]*(Jackpot+1)
    nxtPi[0] = pi[0] + pDn*pi[1]
    nxtPi[1] = pi[2]*pDn

    nxtPi[Jackpot] = pi[Jackpot]+pi[Jackpot-1]*pUp
    nxtPi[Jackpot-1] = pi[Jackpot-2]*pUp

    for i in xrange(2,Jackpot-1):
        nxtPi[i] = pi[i-1]*pUp+pi[i+1]*pDn

    return nxtPi


def init(stake):
    nxtPi = [0]*(Jackpot+1)
    nxtPi[stake] = 1
    return nxtPi

def endpoints(step,pi):
    #print(step,pi[0],pi[Jackpot])
    print "Step %d, P(ruin) = %f, P(jackpot) = %f" % (step, pi[0], pi[Jackpot] )

def converged(xpi,pi):
    maxDelta = 0.0
    for idx in xrange(0,len(xpi)):
        delta = math.fabs(xpi[idx]-pi[idx]) 
        if maxDelta < delta:
            maxDelta = delta

        if 5e-8 < maxDelta:
            return False

    return True

pi = init(Stake)

step = 0
while(True):
    xpi = nxtStep(pi)
    endpoints(step,pi)
    step += 1
    if converged(xpi,pi):
        break
    pi = xpi

