#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script pour se connecteur sur serveur des données ECMWF et télécharger le modèle demandé
# Plus d'informations au niveau des paramètres sur https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
# On peut également constituer la requête grâce au catalogue http://apps.ecmwf.int/data-catalogues/era5/?class=ea
# Afin de calculer le géopotentiel et la pression, il est nécessaire de récupérer la température, l'humidité spécifique, le logarithme de la pression surfacique et le géopotentiel au sol (IDs : 130, 133, 152, 129)

from ecmwfapi import ECMWFDataServer
 
server = ECMWFDataServer()
server.retrieve({
    "class": "ea", # A ne pas changer
    "dataset": "era5", # A ne pas changer
    "expver": "1", # A ne pas changer
    "stream": "oper", # Atmospheric model
    "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137",
    "type": "an", # Type hindcast (analysis) et non forecast
    "levtype": "ml", # Model levels (défini au https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels)
    "param": "129/130/131/132/133/135/152/203", # ID des paramètres, à regarder par exemple dans la Table 12 du https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
    "date": "2016-05-27/to/2016-05-28", # Plage des dates
    # time : default 1hr
    "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
    "step": "0", # Pour hindcast (analysis) doit être à 0 (pour forecast, peut changer)
    "grid": "1.0/1.0", # Default @1/1, sinon interp pour moins (see https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step)
    "area": "-36/-80/-37/-62", # La zone à récupérer : latmin/longmin/latmax/longmax
    "format": "netcdf", # Format des données, netcdf est facilement lu par Matlab
    "target": "era5.nc" # Le nom du fichier obtenu
 })
