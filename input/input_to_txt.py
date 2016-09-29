#!/usr/bin/python

import scipy.io as sio
import scipy;
import numpy as np
import sys
import math
import os.path

if (len(sys.argv) < 2):
	print "Usage: to_txt.py <input_graph1.mat> <input_graph2.mat> ..."
	exit(1)

for fidx in range(1, len(sys.argv)):
	f = sys.argv[fidx]
	if f[-3:] == "txt":
		continue
	if f[-3:] != "mat":
		print "format of file %s is not supported." % (f)
		continue
	if os.path.isfile(f + '.txt'):
		continue

	print "transforming %s to txt file..." % (f)
	m = sio.loadmat(f)
	data = m['Problem'][0][0]
	idx = 0;
	for i in range(len(data)):
		#print ("%d %s" % (i, type(data[i])))
		if (type(data[i]) is scipy.sparse.csc.csc_matrix):
			idx = i;
			break;
	dense_m = data[idx].todense()
	outfile = f + '.txt'
	out = open(outfile, 'w')
	out.write('#Converted from %s\n' % (f))
	out.write('#Nodes: %d\n' % (dense_m.shape[0]))
	it = np.nditer(dense_m, flags=['multi_index'])
	while not it.finished:
		if (math.fabs(it[0]) < 0.00001):
			pass
		else:
			out.write("%d %d %f\n" % (it.multi_index[0], it.multi_index[1], it[0]))
		it.iternext()
	out.close()	
