#!/usr/bin/python

import scipy.io as sio
import scipy;
import numpy as np
import sys
import math

if (len(sys.argv) < 2):
	print "Usage: to_txt.py <input_graph1.mat> <input_graph2.mat> ..."
	exit(1)
for i in range(1, len(sys.argv)):
	f = sys.argv[i]
	m = sio.loadmat(f)
	data = m['Problem'][0][0]
	idx = 0;
	for i in range(len(data)):
		print ("%d %s" % (i, type(data[i])))
		if (type(data[i]) is scipy.sparse.csc.csc_matrix):
			idx = i;
			print "I am here"
			break;
	dense_m = data[idx].todense()
	print dense_m
	outfile = sys.argv[1] + '.txt'
	out = open(outfile, 'w')
	out.write('#Converted from %s\n' % (sys.argv[1]))
	out.write('#Nodes: %d  Edges: %d\n' % (dense_m.shape[0], dense_m.shape[1]))
	it = np.nditer(dense_m, flags=['multi_index'])
	while not it.finished:
		if (math.fabs(it[0]) < 0.00001):
			pass
		else:
			out.write("%d %d %f\n" % (it.multi_index[0], it.multi_index[1], it[0]))
		it.iternext()
	out.close()
