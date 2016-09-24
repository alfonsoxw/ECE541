from __future__ import print_function

import pdb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(u'-busy', metavar  = u'probability channel is busy', dest=u'busy', required=True)
parser.add_argument(u'-fail', metavar = u'fail', dest=u'fail', required=True)
parser.add_argument(u'-steps', metavar = u'fail', dest=u'steps', required=True)
args = parser.parse_args()

channelBusy  = float(args.busy)
msgFail      = (1-channelBusy)*float(args.fail)
steps        = int(args.steps)

p01 = 1
p12 = msgFail
p11 = channelBusy
p20 = 1
p13 = 1-(msgFail+channelBusy)
p33 = 1



def update(pi):
    xpi = [0]*4

    xpi[0] = pi[2]
    xpi[1] = pi[0] + pi[1]*p11
    xpi[2] = pi[1]*p12
    xpi[3] = pi[3] + pi[1]*p13
    return xpi

pi     = [0]*4
pi[0]  = 1

mean = 0.0

for idx in xrange(0,steps):
    pi = update(pi)
    print('step',idx,'P(msg having succeeded by now)',pi[3])
    mean += 1.0-pi[3]


print('estimated mean steps to success',mean)


