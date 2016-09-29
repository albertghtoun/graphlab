#!/usr/bin/python
import sys
from subprocess import call
import os
import os.path

if len(sys.argv) != 2:
	print "usage: download.py <URL file>"
	exit()
fl = open(sys.argv[1], 'r')
for l in fl:
	name = l[l.rindex('/')+1:-1] #remove ending \n
	idx = name.index('.') 
	if idx <= 0:
		print "unknown file type: %s" % (name)
		continue
	files = [f for f in os.listdir('.') if os.path.isfile(f)]
	exist = False
	for f in files:
    		if f[:f.index('.')] == name[:idx]:
			exist = True
			break
	if exist == True:
		continue
	print "downloading %s..." % (name)
	print "wget %s" % (l[:-1])
	call(["wget", l[:-1]])
	if name.endswith('.tar.gz'):
		call["tar", "xzf", name]
		call["rm", "-fr", name]
	elif name.endswith('.txt'):
		pass
	elif name.endswith('.mat'):
		pass
fl.close()
