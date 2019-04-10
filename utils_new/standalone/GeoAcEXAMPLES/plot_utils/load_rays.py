# coding=utf-8
# Author:        Léo Martire.
# Description:   Loads some raypaths previously generated by GeoAc.
# Notes:         
#
# Usage:
#   
# with:
#   
# yields:
#   

import numpy as np # mandatory
import os # nice prints
import re # mandatory

def load_rays(filepath, renorm_dist):
  idsxyz=[0,1,2]
  idsattenuation=[3,4]
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
  # remove eventual empty rays
  for r in range(0,len(fulldata)):
    if(len(fulldata[r])==0):
      fulldata.pop(r)
  NRAYS=len(fulldata)
  maxz=0.
  maxd=0.
  print("["+os.path.basename(__file__)+"] Found "+str(NRAYS)+" rays:")
  for r in range(0,NRAYS):
    print("["+os.path.basename(__file__)+"]   Ray "+str(r+1)+': '+str(np.shape(fulldata[r])))
    # make numpy array
    fulldata[r]=np.array(fulldata[r])
    # convert distances to good units
    fulldata[r][:,idsxyz] *= renorm_dist
    # find eventual nans in attenuation
    if(np.any(np.isnan(fulldata[r][:,idsattenuation]))):
      print("["+os.path.basename(__file__)+"]     Ray "+str(r+1)+": found NaNs in attenuation columns. Replacing them with zeros.")
      fulldata[r][:,idsattenuation]=np.nan_to_num(fulldata[r][:,idsattenuation])
    locmaxd=np.max(np.linalg.norm(fulldata[r][:,[0,1]],axis=1))
    if(locmaxd>maxd):
      maxd=locmaxd
    locmaxz=np.max(fulldata[r][:,2])
    if(locmaxz>maxz):
      maxz=locmaxz
  return(fulldata, NRAYS, maxd, maxz)
