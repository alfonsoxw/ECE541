from __future__ import print_function

import pdb
import argparse
import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.stats import poisson


def ctmcPI(Q,alpha):
	norm = 0.0

	# compute normalization constant
	#
	pi = alpha.copy()
	for idx in xrange(0,5):
		pi[idx,0] = -alpha[idx,0]/Q[idx,idx]
		norm += -alpha[idx,0]/Q[idx,idx]

	# normalize each product of DTMC probability times mean time in state
	for idx in xrange(0,5):
		pi[idx,0] /= norm

	return pi

def convergence(L):
	myL = sorted(L,key=lambda v: np.absolute(v),reverse=True)
	return np.absolute(myL[1]/myL[0])


# main body of code.  Get the arguments
#
parser = argparse.ArgumentParser()
parser.add_argument(u'-request', metavar = u'request rate', dest=u'request', required=True)
parser.add_argument(u'-retry', metavar   = u'communication retry rate', dest=u'retry', required=True)
parser.add_argument(u'-network', metavar = u'effective network rate', dest=u'network', required=True)
parser.add_argument(u'-service', metavar = u'effective service rate', dest=u'service', required=True)
parser.add_argument(u'-networkFail', metavar = u'network failure probability', dest=u'networkFail', required=True)
parser.add_argument(u'-serviceFail', metavar = u'service failure probability', dest=u'serviceFail', required=True)
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

error    = float(args.error) 

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


Pt = P.T

# create embedded Markov chain
#
eP = Q.copy()
for rdx in xrange(0,5):
	diag = -Q[rdx,rdx]
	for cdx in xrange(0,5):
		eP[rdx,cdx] /= diag
	eP[rdx,rdx] = 0.0
	


ePt = eP.T
 
x   = np.matrix([0.2, 0.2, 0.2, 0.2, 0.2]).T
ex  = np.matrix([0.2, 0.2, 0.2, 0.2, 0.2]).T

k=1
dnorm = 2*error
enorm = dnorm

while error < dnorm and error < enorm:

	# apply power law to I+Q/lambda matrix

	nx = Pt*x
	dx = x-nx
	dnorm = np.linalg.norm(dx)
	x = nx

	# apply power law to embedded DTMC probability vector
	#
	enx = ePt*ex
	dx  = ex-enx
	enorm = np.linalg.norm(dx)
	ex = enx

	k += 1
	if not k%10:
		print(k, dnorm, enorm)


# direct solve, to check power method answer
#
tQ = Q.copy().T

# force row 0 to be all ones
#
for idx in xrange(0,5):
	tQ[0,idx] = 1.0

# unit column vector
#
e0 = np.array([1.0,0,0,0,0]).T

# direct solver in numpy uses LU decomposition
#
direct_pi = np.linalg.solve(tQ,e0)

ex = ctmcPI(Q,ex)

print('direct ', direct_pi[0],direct_pi[1],direct_pi[2],direct_pi[3],direct_pi[4])
print('power p', x[0,0],x[1,0],x[2,0],x[3,0],x[4,0])
print('power e', ex[0,0],ex[1,0],ex[2,0],ex[3,0],ex[4,0])

eigP = np.linalg.eig(Pt)
eigE = np.linalg.eig(ePt)

# eigP[0] is list of eigenvalues.
#
convergenceRateP = convergence(eigP[0])
convergenceRateE = convergence(eigE[0])

print('scaled Q convergence',convergenceRateP,'embedded DTMC convergence', convergenceRateE)



