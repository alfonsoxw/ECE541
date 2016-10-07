from __future__ import print_function

import pdb
import argparse
import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.stats import poisson

# main body of code.  Get the arguments
#
parser = argparse.ArgumentParser()
parser.add_argument(u'-request', metavar = u'request rate', dest=u'request', required=True)
parser.add_argument(u'-retry', metavar   = u'communication retry rate', dest=u'retry', required=True)
parser.add_argument(u'-network', metavar = u'effective network rate', dest=u'network', required=True)
parser.add_argument(u'-service', metavar = u'effective service rate', dest=u'service', required=True)
parser.add_argument(u'-networkFail', metavar = u'network failure probability', dest=u'networkFail', required=True)
parser.add_argument(u'-serviceFail', metavar = u'service failure probability', dest=u'serviceFail', required=True)
parser.add_argument(u'-time', metavar = u'time t in transient solution', dest=u'time', required=True)
parser.add_argument(u'-error', metavar = u'error bound', dest=u'error', required=True)
args = parser.parse_args()

# compute transition rates from parameters
#
lambdaF = float(args.networkFail)*float(args.network)
lambdaN = float(args.network) - lambdaF

# compute transition rates from parameters
#
lambdaSF = float(args.serviceFail)*float(args.service)
lambdaS  = float(args.service) - lambdaSF

lambdaT  = float(args.request)
lambdaR  = float(args.retry)

time     = float(args.time)  if args.time else None
error    = float(args.error) if args.error else None

# create array for rate matrix
#
zeros = [0.0]*5
zerosList = []
for idx in xrange(0,5):
    zerosList.append(zeros)

Q = np.matrix(zerosList)

# fill in rate structure from diagram
#
Q[0,1] = lambdaT
Q[0,0] = -lambdaT

Q[1,2] = lambdaN
Q[1,3] = lambdaF
Q[1,1] = -(lambdaN+lambdaF)

Q[2,4] = lambdaS
Q[2,3] = lambdaSF
Q[2,2] = -(lambdaS+lambdaSF)

Q[3,1] = lambdaR
Q[3,3] = -lambdaR

Q[4,0] = lambdaN
Q[4,3] = lambdaF
Q[4,4] = -(lambdaN+lambdaF)


# direct solve, to check uniformation answer
#
tQ = Q.copy().T

# force row 0 to be all ones
#
for idx in xrange(0,5):
	tQ[0,idx] = 1.0

# unit column vector
#
e0 = np.array([1.0,0,0,0,0]).T


direct_pi = np.linalg.solve(tQ,e0)

# find largest in magnitude diagonal element
#
unifRate = 0.0
for idx in xrange(0,5):
	if unifRate < -Q[idx,idx]:
		unifRate = -Q[idx,idx]


# transform Q to I + Q/unifRate 
#
P = Q.copy()

for idx in xrange(0,5):
	for jdx in xrange(0,5):
		# scale every element by 1/unifRate
		P[idx,jdx] /= unifRate

		# add 1 on diagonal
		if idx==jdx:
			P[idx,idx] += 1.0


# accumProb will accumulate what amounts to P( Poisson(unifRate) <= k )
#
accumProb = math.exp(-unifRate*time)

# DTMC step
k = 1

# accumulate pi(t) in array
#pi  = [ math.exp(-unifRate*time), 0.0, 0.0, 0.0, 0.0]
pi  = np.array([math.exp(-unifRate*time), 0.0, 0.0, 0.0, 0.0])

#   x will be probability occupancy vector for DTMC, a numpy
#   'array' so that we can call numpy to do vector-matrix multiplication
#        for x(k) = x(k-1) P
#

x   = np.array([1.0, 0.0, 0.0, 0.0, 0.0])

# while the error bound exceeds the error threshold
#  keep iterating

while( error < 1.0-accumProb):

	# advance DTMC probabilities via x(k) = x(k-1) P
	#
	x = x.dot(P)

	# get probability of k DTMC transitions in uniformized chain
	#
	prK = poisson.pmf(k,unifRate*time)

	# record this accumulation for error bound
	accumProb += prK

	# add to pi the contribution made by k transitions of the DTMC
	#
	pi = pi + prK*x

	#for idx in xrange(0,5):
	#	pi[idx] = pi[idx] + prK * x[0,idx]

	# one more step
	#
	k += 1
	if not k%20:
		print(k,1.0-accumProb)

print('transient ', pi[0,0],pi[0,1],pi[0,2],pi[0,3],pi[0,4])
print('stationary', direct_pi[0],direct_pi[1],direct_pi[2],direct_pi[3],direct_pi[4])

