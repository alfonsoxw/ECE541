from __future__ import print_function
import pdb
import math
import StringIO


N = 10 
Start = 4

nth = 1.0/N


def nxtStep(pi):
    nxtPi = [0]*(N+1)
    nxtPi[0] = pi[1]*nth
    nxtPi[N] = pi[N-1]*nth
   
    for i in xrange(1,N):
        pUp = (N-i+1)*nth   # chance of selecting a ball from 'other' urn in state i-1
        pDn = (i+1)*nth     # chance of selecting ball from tracked urn in state i+1

        nxtPi[i] = pi[i-1]*pUp + pi[i+1]*pDn

    return nxtPi


def init(cnt):
    nxtPi = [0]*(N+1)
    nxtPi[cnt] = 1
    return nxtPi

def endpoints(step,pi):
    #print(step,pi[0],pi[Exit])
    out = StringIO.StringIO()
    for v in pi:
        out.write( '%f ' % v )
    print('Step',step,out.getvalue()) 

def converged(xpi,pi):
    maxDelta = 0.0
    for idx in xrange(0,len(xpi)):
        delta = math.fabs(xpi[idx]-pi[idx]) 
        if maxDelta < delta:
            maxDelta = delta

        if 5e-8 < maxDelta:
            return False

    return True

pi = init(Start)
oldpi = [0]*(N+1)
step = 0
while(True):
    xpi = nxtStep(pi)
    endpoints(step,pi)
    step += 1
    if converged(xpi,oldpi):
        break

    oldpi = pi
    pi = xpi

