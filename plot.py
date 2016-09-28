"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt

def load_data(datafile):
	data = {}
	f = open(datafile, 'r')
	header = True
	for line in f:
		stringlist = line.split()
		if header == True:
			for i in range(len(stringlist)):
				data[stringlist[i]]=[]
			header=False
			continue

n_groups = 10
data = load_data()
data = (20, 35, 30, 35, 27, 20, 30, 40, 49, 10)

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.08
opacity = 0.4

for i in range(n_groups):
	rects2 = plt.bar(index + i*bar_width, data, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Women')

plt.xlabel('Gather Percent(%)')
plt.ylabel('Time(s)')
plt.title('Running time under different gather percent')
plt.xticks(index + bar_width, ('10', '20', '30', '40', '50', '60', '70', '80', '90', '100'))
plt.legend()

plt.tight_layout()
plt.show()
