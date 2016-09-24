from __future__ import print_function

import pdb
import argparse
import numpy as np
from scipy.stats import poisson
from math import fabs
import matplotlib.pyplot as plt
from discreteMarkovChain import markovChain

# Markov chain model of inventory model, known in the literature as (s,S) inventory model.
#
# Every day in a store starts with some K items.  A random number of items are purchased
# during the day, the distribution of the number _attempted_ for purchase is i.i.d every day, in this
# solution, Poisson.
#
# The maximum inventory size is N.   If the inventory drops to or below 'restock',
# then at the end of the day the inventory will be refiled to its maximum level, N.
# If the demand during a day does not drop the level to or below restock, it remains the same to the next day.
#
# The state of the chain is the inventory level at the beginning of the day...note that it will reflect
# any restocking that occured at the end of the prior day.
#


# probability of making N to N transition.  N is special because you can get to 
# N in one day two different ways, one when you start with N and there is no demand,
# the other when you start in N and have enough demand to restock to N
#
def probNtoN():

    # chance of no demand
    #
    prob  =  poisson.pmf(0,dr)                       # note that poisson.pmf is probability mass function of 
                                                     # Poisson, with rate 'dr'

    # add chance of getting a large enough demand to restock
    #
    prob +=  (1.0 - poisson.cdf(N-(resupply+1),dr))  # note that poisson.cdf is cumulative distribution function

    return prob

# probability of making J to N transition, with J<N.  
# Only way this happens is with demand large enough to
# drop inventory level to or below resupply level 
#
def probJtoN(J):
    return (1.0-poisson.cdf(J-(resupply+1),dr))


# probability of making J to K transition, with K <=J < N.
# Only way this happens is if exactly J-K items are purchased.
# Function only called on K larger than resupply
#
def probJtoK(J,K):
    return poisson.pmf(J-K,dr)

# we don't use states for inventory levels resupply or less,
# so provide a mapping from 'Inventory' coordinates to 'array' coordinates
#
def I2a(i):
    return i-(resupply+1) 


# runtime parameters are restocking level, rate parameter for Poisson distribution, and inventory capacity
#
parser = argparse.ArgumentParser()
parser.add_argument(u'-s', metavar  = u'restocking level', dest=u's', required=True)
parser.add_argument(u'-dr', metavar = u'Demand rate', dest=u'dr', required=True)
parser.add_argument(u'-N', metavar = u'Capacity', dest=u'N', required=True)
args = parser.parse_args()

resupply  = int(args.s)
N         = int(args.N)
dr        = float(args.dr)



# Create the transition matrix.  Problem is small enough that there's no need to
# be fancy, just initialize a full matrix with zeros
#
States    = N-resupply
zeros     = [0.0]*States     # initialize an array with each element having value 0.0

zerosList = []               # a matrix is initialized by a list of lists

for n in xrange(0,States):   # make a list with States elements, each component being a list with States zeros
   zerosList.append(zeros)

# Need to use numpy matrix when we do a solve on the transition matrix
P 	  = np.matrix(zerosList)


# flesh out the transition matrix
for J in xrange(resupply+1, N+1):

	# probability of transitioning up to full inventory
	#
	P[ I2a(J), I2a(N) ] = probJtoN( J )

	# probability of transitioning to a lower level without restocking
	#
	for K in xrange(resupply+1, J+1):
	   P[ I2a(J), I2a(K) ] = probJtoK( J, K )

# overwrites same entry computed in above loop, the Nth-Nth element is special
#
P[ I2a(N) ,I2a(N) ] = probNtoN()

# check the transition matrix for double stochasticity
#  First the rows, each has to sum to 1
#
for i in xrange(resupply+1,N+1):
	sum = 0.0
	for j in xrange(resupply+1,N+1):
		sum += P[I2a(i), I2a(j)]
	if 1e-9 < fabs(1.0-sum):
		print('suspect transition matrix row')
		print(P)
		exit(0)	

# use the discrete state markov solver to get the steady state probabilities
#
mc = markovChain(P)
mc.computePi('linear')

# mc.pi is the solution vector, print it out
print (mc.pi)

