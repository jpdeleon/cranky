# general python files
from os import listdir
from os.path import join
from glob import glob
import numpy as np
import getpass
from run_pipeline import run


inputpath='new_example_input/'
# get all files in the folder, just grab the EPIC number of the filename
fitfiles = glob(join(inputpath,'*.fits'))
starnames = np.unique([f.split('/')[-1][4:13] for f in fitfiles])

#inputpath='/home/'+getpass.getuser()+'/data/transit/c14/'
#fitfiles = glob(join(inputpath,'*.fits'))
#starnames = np.unique([f.split('/')[-1][4:13] for f in fitfiles])

outputpath='new_example_output/'
#print(fitfiles)

i = 0

exc_list = []
while i < len(starnames):
  print('Now running stars, number ')
  print(str(i))

  try:
    run(starnames[i],outputpath=outputpath,inputpath=inputpath,makelightcurve=True, showfig=True)#,campaign=14)

  except Exception as inst:
    print inst
    exc_list.append(inst)
  i = i + 1

if exc_list:
  print('Module failed for some stars:')
  print(exc_list)

print('Done...')
