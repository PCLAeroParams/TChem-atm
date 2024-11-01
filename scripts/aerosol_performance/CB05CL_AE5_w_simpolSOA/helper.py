import glob, os, json
import pandas as pd
import numpy as np
def create_dataFrame_cpu(flist):
  df_list = []
  for file in flist:
    fname = os.path.basename(file)
    # fileNames.append(fname)
    colName = np.array(fname.replace('.json', '').split('_'))
    nBatch = colName[3]
    nThread = colName[6]
    nParticles =colName[10]
    colName = '_' + nBatch + 'b_' + nThread + 't' + nParticles +'p'
    data = pd.read_json(file)
    # data.rename(columns={'Aerosol RHSs': 'RHSwalltime'}, inplace=True)
    # data.rename(columns={'Aerosol Numerical Jacobian': 'Jacwalltime'}, inplace=True)
    data.rename(lambda x: x.replace('wall_time_', ''), axis='index', inplace=True)
    data.rename({'number_of_samples': 'nSamples', 'number_of_time_iters': 'nIter',
                 'iter_0': 'walltime'},
                  axis='index', inplace=True)
    data.loc['nBatch'] = int(nBatch)
    data.loc['nThread'] = int(nThread)
    data.loc['nParticles'] = int(nParticles)
    tmp = data.add_suffix(colName, axis=1)
    df_list.append(tmp)
    df = pd.concat(df_list, axis=1)
    df.sort_values(by=['nBatch', 'nThread','nParticles'], axis=1, inplace=True)
    df = df.transpose()
  return df


def create_dataFrame_gpu(flist):
  df_list = []
  for file in flist:
    fname = os.path.basename(file)
    # fileNames.append(fname)
    try:
      data = pd.read_json(file)
      colName = np.array(fname.replace('.json', '').split('_'))
      # print(colName)
      nBatch = colName[3]
      # note: we have this set to 1 currently, but keeping in case of change
      vector_size = colName[5]
      team_size = colName[8]
      nParticles =colName[12]
      # colName = '_' + nParticles +'p'
      colName = '_' + nBatch+'t_'+ team_size  +'v_'+ vector_size  + 'b_'  + nParticles +'p'

      data.rename(columns={'Atmospheric Chemistry E3SM': 'walltime'}, inplace=True)
      data.rename(lambda x: x.replace('wall_time_', ''), axis='index', inplace=True)
      data.rename({'number_of_samples': 'nSamples', 'number_of_time_iters': 'nIter',
                 'iter_0': 'walltime'},
                  axis='index', inplace=True)
      data.loc['nBatch'] = int(nBatch)
      data.loc['vector_size'] = int(vector_size)
      data.loc['team_size'] = int(team_size)
      data.loc['nParticles'] = int(nParticles)
      tmp = data.add_suffix(colName, axis=1)
      df_list.append(tmp)
    except:
      print(file + "has issues.")
  df = pd.concat(df_list, axis=1)
  df.sort_values(by=['nBatch', 'vector_size','team_size','nParticles'], axis=1, inplace=True)
  df = df.transpose()
  return df
