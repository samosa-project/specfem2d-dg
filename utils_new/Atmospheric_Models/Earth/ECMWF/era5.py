#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:        Guerman Poler, Léo Martire.
# Description:   Calls ECMWF API with requested parameters.
# Last modified: See file metadata.
# Usage:         Call this script with 5 arguments (startDate, endDate, times, areaStr, outputFileName).
# Notes:         More information on parameters at https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation.
#                Building the request can also be done using catalog at http://apps.ecmwf.int/data-catalogues/era5/?class=ea.
#                In order to compute geopotential and pressure at all levels, one needs to request temperature, specific humidity, logarithm of surface pressure, and ground geopotential (IDs 130, 133, 152, and 129 as defined in Table 12 at https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation).

from ecmwfapi import ECMWFDataServer
import sys

#print 'Argument list:', str(sys.argv) # DEBUG

script_name=((sys.argv[0].split('/'))[1].split('.'))[0]

if(len(sys.argv)<6):
  print('['+script_name+', ERROR] Not enough input arguments. I need at least 5 (startDate, endDate, times, areaStr, outputFileName).')
  exit(-1)

# Assign variables from arguments.
startDate=sys.argv[1]
endDate=sys.argv[2]
time=sys.argv[3]
area=sys.argv[4]
outputFile=sys.argv[5]

dateStr=(startDate+"/to/"+endDate)

print('['+script_name+'] Calling API with parameters:')
print('['+script_name+']   date:   '+dateStr)
print('['+script_name+']   time:   '+time)
print('['+script_name+']   area:   '+area)
print('['+script_name+']   target: '+outputFile)

exit()

server = ECMWFDataServer()

server.retrieve({
    # Do not change.
    "class":    "ea",
    # Do not change.
    "dataset":  "era5",
    # Do not change.
    "expver":   "1",
    # Do not change: asks for an atmospheric model.
    "stream":   "oper",
    # Do not change: requests all levels (1-137).
    "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137",
    # Do not change: asks for hindcast (analysis) and not forecast.
    "type":     "an",
    # Do not change: model levels (as defined at https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels).
    "levtype":  "ml",
    # Do not change: parameters' IDs (as defined in Table 12 at https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation).
    "param":    "129/130/131/132/133/135/152/203",
    
    # Dates. Format YYYY-MM-DD/to/YYYY-MM-DD.
    #"date":     "2016-05-27/to/2016-05-27",
    "date":     dateStr,
    
    # Times for each day (same for each day). Format HH:MM:SS[/HH:MM:SS[/HH:MM:SS[/...]]]. By default, step is 1 hour. Asking finer may result in interpolated values.
    #"time":     "00:00:00",
    "time":     time,
    
    # Do not change: for hindcasts (analysis) must be 0 (for forecasts, may change).
    "step":     "0",
    # Do not change: spatial grid step. By default, step is 1°/1° (1.0/1.0). Finer may result in interpolated results (see https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step).
    "grid":     "1.0/1.0",
    
    # Area. Format latmin/longmin/latmax/longmax.
    #"area":     "-36/-80/-37/-62",
    "area":     area,
    
    # Do not change: data format.
    "format":   "netcdf",
    
    # Name of output file.
    #"target":   "era5.nc"
    "target":   outputFile
 })
