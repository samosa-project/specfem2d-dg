import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from mpl_toolkits.mplot3d import Axes3D  # for 3d plot
import matplotlib.collections as mcoll # for gradient-colored lines
import matplotlib.path as mpath # for gradient-colored lines

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
rc('font', **font)
rc('text', usetex=True)

###############################
def load_rays(filepath, renorm_dist):
  fulldata=[]
  ray=1
  fulldata.append([]) # add new list for new ray
  with open(filepath, 'rt') as f:
    line=f.readline()
    while line:
      p=re.compile('^ *#')
      #print(line)
      #print(p.findall(line))
      if(len(p.findall(line))==0):
        #print('non-comment line found')
        p=re.compile('^ *-?[0-9]+')
        #print(p.findall(line))
        if(len(p.findall(line))>0):
          #print('number line found')
          #print(line)
          fulldata[ray-1].append(np.fromstring(line, sep=' '))
        else:
          #print('space line found')
          ray+=1
          fulldata.append([]) # add new list for new ray
      line=f.readline()

  for r in range(0,len(fulldata)):
    #print(len(fulldata[r]))
    if(len(fulldata[r])==0):
      fulldata.pop(r) # remove eventual empty rays
    else:
      # make numpy array
      fulldata[r]=np.array(fulldata[r])
      # convert distances to good units
      fulldata[r][:,0:3] *= renorm_dist
      # find eventual nans in attenuation
      

  NRAYS=len(fulldata)
  print(str(NRAYS)+" rays found:")
  maxz=0.
  maxd=0.
  for r in range(0,NRAYS):
    print('  ray '+str(r+1)+': '+str(np.shape(fulldata[r])))
    locmaxz=np.max(fulldata[r][:,2])
    locmaxd=np.max(np.linalg.norm(fulldata[r][:,[0,1]],axis=1))
    if(locmaxd>maxd):
      maxd=locmaxd
    if(locmaxz>maxz):
      maxz=locmaxz
  #print(fulldata)
  #print(len(fulldata))
  #print(np.shape(fulldata[1])[1])
  #print(np.shape(fulldata[0][0]))
  return(fulldata, NRAYS, maxd, maxz)
###############################
def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()
    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)
    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])
    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
###############################
def make_segments(x, y):
  """
  Create list of line segments from x and y coordinates, in the correct format
  for LineCollection: an array of the form numlines x (points per line) x 2 (x
  and y) array
  """
  points = np.array([x, y]).T.reshape(-1, 1, 2)
  segments = np.concatenate([points[:-1], points[1:]], axis=1)
  return segments
def colorline(
  x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
      linewidth=3, alpha=1.0):
  """
  http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
  http://matplotlib.org/examples/pylab_examples/multicolored_line.html
  Plot a colored line with coordinates x and y
  Optionally specify colors in the array z
  Optionally specify a colormap, a norm function and a line width
  """
  # Default colors equally spaced on [0,1]:
  if z is None:
      z = np.linspace(0.0, 1.0, len(x))
  # Special case if a single number:
  if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
      z = np.array([z])
  z = np.asarray(z)
  segments = make_segments(x, y)
  lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                            linewidth=linewidth, alpha=alpha)
  ax = plt.gca()
  ax.add_collection(lc)
  return lc
###############################

parser = argparse.ArgumentParser(description='Plots some raypaths.')

required=parser.add_argument_group('required arguments')
required.add_argument("-f", "--file", type=str, help="file containing raypaths", required=True)

parser.add_argument("-Z", "--zmax", type=float, default=float('Inf'), help="z_max for plot")
parser.add_argument("-d", "--dmin", type=float, default=0., help="d_min for plot")
parser.add_argument("-D", "--dmax", type=float, default=float('Inf'), help="d_max for plot")
parser.add_argument("-p", "--projType", type=str, default='2d', help="Projection type.", choices=['2d', '3d'])
parser.add_argument("-a", "--maxAttenuation", type=float, default=-50., help="Maximum attenuation (dB).")

args = parser.parse_args()

unit_import_intermsofmeters=1000 # outputfile specifies km
unit_plot_intermsofmeters=1 # we want to plot in m

# load
(fulldata, NRAYS, datmaxr, datmaxz) = load_rays(args.file, unit_import_intermsofmeters / unit_plot_intermsofmeters)
exit()

#print(np.sum(fulldata[0][0:3,0:2],axis=1)); exit()

maxz=min(args.zmax,datmaxz)
minr=max(args.dmin,0.)
maxr=min(args.dmax,datmaxr)
plottype=args.projType
maxattenuationdb=args.maxAttenuation

idx=0
idy=1
idz=2
idattgeo=3
idattatm=4
idt=5
idd=np.shape(fulldata[0])[1]

# plot choices in terms of unit of plot
#plottype='2D'
#plottype='3D'

fig = plt.figure()
if(plottype=='3d'):
  ax = fig.add_subplot(111, projection='3d')
  ptype=3
elif(plottype=='2d'):
  ax = fig.add_subplot(111, aspect='equal')
  ptype=2
else:
  print('error'); exit(-1)

for r in range(0,NRAYS):
#for r in [1]:
  ray=fulldata[r]
  d=np.linalg.norm(ray[:,[idx,idy]],axis=1)
  sel=np.logical_and(d>minr,d<maxr)
  ncolz=np.shape(ray)[1]
  nlines_sel=np.sum(sel)
  selectedray=np.zeros((nlines_sel,ncolz))
  for c in range(0,ncolz):
    selectedray[:,c]=np.extract(sel,ray[:,c])
  seld=np.reshape(np.extract(sel,d),(nlines_sel,1))
  #print(np.shape(selectedray))
  #print(np.shape(seld))
  #print(selectedray)
  #print(seld)
  selectedray=np.append(selectedray, seld,axis=1)
  #print(np.shape(selectedray))
  #print(selectedray)
  #print(selectedray[:,idx])
  colour = np.extract(sel,ray[:,idt])/520.
  colour = np.extract(sel,np.sum(ray[:,[idattgeo,idattatm]],axis=1))/maxattenuationdb
  #print('minmax col'+str(min(colour))+' '+str(max(colour)))
  #print(colour)
  
  if(ptype==3):
    #ax.scatter(ray[:,idx],ray[:,idy],ray[:,idz])
    ax.scatter(selectedray[:,idx],selectedray[:,idy],selectedray[:,idz])
  elif(ptype==2):
    #ax.plot(selectedray[:,idd],selectedray[:,idz])
    colorline(selectedray[:,idd], selectedray[:,idz], colour, cmap=plt.get_cmap('gist_heat'), linewidth=2)
  else:
    print('error'); exit(-1)

plt.axis('scaled')
#set_aspect_equal_3d(ax)

if(ptype==3):
  ax.set_zlim(0, maxz)
elif(ptype==2):
  ax.set_ylim(0, maxz)
else:
  print('error'); exit(-1)

#plt.colorbar()
plt.grid(True)
plt.show()
