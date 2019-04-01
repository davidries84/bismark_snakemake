"""
Histogram of samtools depth output


    * Setting the number of data bins
    * The ``normed`` flag, which normalizes bin heights so that the integral of
      the histogram is 1. The resulting histogram is a probability density.
    * Setting the face color of the bars
    * Setting the opacity (alpha value).

"""
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

# load data
f = open( sys.argv[1], 'r')
x = []
for line in f.readlines():
	x.append(int(line.strip().split()[2]))

f.close()

x = np.asarray(x)

mean = np.mean(x)
median = np.median(x)

print mean
print median

num_bins = 200
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, range=(0,200))#, normed=1, facecolor='green', alpha=0.5)

# plot mean
plt.plot((mean, mean),(0,max(n)), 'k-')
plt.plot((median, median), (0, max(n)), 'b-')
plt.annotate('mean '+str(mean), xy=(mean, max(n)))
plt.annotate('median '+str(median), xy=(median, max(n)))

plt.xlabel('Coverage')
plt.ylabel('# of Bases')
plt.title(r'Histogram of coverage')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()
