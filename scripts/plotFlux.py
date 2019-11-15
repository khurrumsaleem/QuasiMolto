# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import fnmatch
import os
import csv
import matplotlib.colors as colors

filePrefix='scalar-flux-group-'

numGroups=len(fnmatch.filter(os.listdir("."), filePrefix+ '*'))

with open('r-mesh.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter=',')
  for row in readCSV:
    R = np.array(list(map(float, row)))

with open('z-mesh.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter=',')
  for row in readCSV:
    Z = np.array(list(map(float, row)))

fig = plt.figure()
fig, axes = plt.subplots(nrows=1, ncols=numGroups)
Z, R = np.meshgrid(Z, R)
flux=np.zeros([len(Z),len(R),numGroups])

for iGroup in np.arange(0,numGroups):
  fileName=filePrefix+str(iGroup)+'.csv'
  with open(fileName) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for count,row in enumerate(readCSV):
      flux[count,:,iGroup] = np.array(list(map(float, row)))

minVal=np.min(flux)
maxVal=np.max(flux)

# Plot the surface.
if (numGroups > 1):
  for count,ax in enumerate(axes.flat):
    colorGrid=ax.pcolormesh(Z,R,flux[:,:,count],cmap=cm.viridis,
                          norm=(colors.LogNorm(vmin=minVal,vmax=maxVal)))
#  colorGrid=ax.pcolormesh(Z,R,flux[:,:,count],cmap=cm.viridis,
#                          vmin=minVal,vmax=maxVal);
    ax.set_title('Group '+str(count))
  fig.colorbar(colorGrid, ax=axes.ravel().tolist(),extend='max')
else:
  colorGrid=axes.pcolormesh(Z,R,flux[:,:,0],cmap=cm.viridis,
                          norm=(colors.LogNorm(vmin=minVal,vmax=maxVal)))

  axes.set_title('Group '+str(0))

  fig.colorbar(colorGrid, ax=axes,extend='max')

fig.suptitle('Scalar Flux',fontsize=16)
plt.savefig("flux.png")
