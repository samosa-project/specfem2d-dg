#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def qualityControlOneSource(_source):
  source = _source.copy()
  
  # Check required keys.
  req_keys = ['xs', 'zs', 'time_function_type', 'factor']
  def_keys = {}
  for rk in req_keys:
    if(not rk in source):
      sys.exit('[%s] Please provide \'%s\' for current source.' % (sys._getframe().f_code.co_name, rk))
  
  # Validate source_type (against Mxx/Mxz/Mzz).
  allMxxMxzMzz = ('Mxx' in source and 'Mxz' in source and 'Mzz' in source)
  anyMxxMxzMzz = ('Mxx' in source or 'Mxz' in source or 'Mzz' in source)
  if('source_type' in source):
    # Check.
    if(source['source_type'] == 2 and not allMxxMxzMzz):
      sys.exit('source_type ==2, but not all of (Mxx, Mxz, Mzz) were found. Something is wrong.')
    if(source['source_type'] == 1 and anyMxxMxzMzz):
      sys.exit('source_type == 1, and some of (Mxx, Mxz, Mzz) were found. Please do not provide them when source_type == 1.')
    else:
      for mmm in ['Mxx', 'Mxz', 'Mzz']:
        source[mmm] = 0.
  else:
    # Try to fix it.
    if(anyMxxMxzMzz):
      if(not allMxxMxzMzz):
        sys.exit('One of Mxx, Mxz, Mzz found, but not all of them. Something is wrong.')
      else:
        source['source_type'] = 2
    else:
      source['source_type'] = 1
  
  # Validate time_function_type.
  if(source['time_function_type']==8):
    if(not 'name_of_source_file' in source):
      sys.exit('Please provide name_of_source_file if time_function_type == 8.')
  else:
    if('name_of_source_file' in source):
      sys.exit('time_function_type != 8, but name_of_source_file was found. Please do not provide it when time_function_type != 8.')
    else:
      source['name_of_source_file'] = '\'NONE\''
  
  # Just set burst_band_width.
  source['burst_band_width'] = 200.
  
  # Validate time_function_type.
  if(source['time_function_type'] in [1,2,3]):
    if(not 'f0' in source):
      sys.exit('Please provide \'f0\' if time_function_type==%d.' % (source['time_function_type']))
  else:
    if('f0' in source):
      sys.exit('Please do not provide \'f0\' if time_function_type==%d.' % (source['time_function_type']))
    else:
      source['f0'] = 1
  
  # Just set tshift and anglesource.
  source['tshift'] = 0.
  source['anglesource'] = 0.
  
  return(source)

def writeOneSource(file_w, source):
  keys = ['xs', 'zs',
          'source_type', 'time_function_type', 'name_of_source_file', 'burst_band_width',
          'f0', 'tshift', 'anglesource',
          'Mxx', 'Mxz', 'Mzz',
          'factor']
  file_w.write('source_surf         = .false.\n')
  for k in keys:
    file_w.write('%-19s = ' % (k))
    if(k=='source_type' or k=='time_function_type'):
      file_w.write('%d' % (source[k]))
    elif(type(source[k])==str):
      file_w.write('%s' % (source[k]))
    else:
      file_w.write('%16.9e' % (source[k]))
    if(k+'_COMMENT' in source):
      file_w.write(' # %s' % (source[k+'_COMMENT']))
    file_w.write('\n')

def roundFieldsToZero(_source):
  source = _source.copy()
  keys = ['Mxx', 'Mxz', 'Mzz',
          'factor']
  for k in keys:
    if(abs(source[k])<1e-16):
      source[k] = 0.0
  return(source)

def writeSources(filepath, sourcesList):
  with open(filepath, 'w') as file_w:
    for isource, source in enumerate(sourcesList):
      source = qualityControlOneSource(source)
      source = roundFieldsToZero(source)
      file_w.write('# SOURCE %04d ########################\n' % (isource+1))
      writeOneSource(file_w, source)
      file_w.write('######################################\n')
      file_w.write('\n')