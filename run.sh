#!/bin/bash

for file in `ls input/*.txt`
do
	./debug/apps/pagerank/pagerank $file
done
