#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script pour se connecteur sur serveur des données ECMWF et télécharger le modèle demandé
# Plus d'informations au niveau des paramètres sur https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
# On peut également constituer la requête grâce au catalogue http://apps.ecmwf.int/data-catalogues/era5/?class=ea
# Afin de calculer le géopotentiel et la pression, il est nécessaire de récupérer la température, l'humidité spécifique, le logarithme de la pression surfacique et le géopotentiel au sol (IDs : 130, 133, 152, 129)

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
    # TODO: explain.
    "levelist": "1",
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
