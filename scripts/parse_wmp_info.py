#!/usr/bin/env python3

import os
import glob
import numpy as np

import openmc.data

wmp_dir = "../WMP_Library" # WMP library PATH
wmp_files = glob.glob(os.path.join(wmp_dir, "*.h5"))

nuclides = []
headers = [
          'Nuclide',
          'WMP File',
          'Energy Range',
          '# Poles',
          '# Windows',
          'CF Order',
          'By VF'
         ]

for i, wmp_library in enumerate(wmp_files):
  result = []
  wmp = openmc.data.WindowedMultipole.from_hdf5(wmp_library)

  wmp_name = os.path.basename(wmp_library)
  atomic_number = int(wmp_name[0:3])
  mass_number = int(wmp_name[3:6])
  nuc_name = openmc.data.data.ATOMIC_SYMBOL[atomic_number] + str(mass_number)
  if wmp_name[6:-3] is not '':
    nuc_name += '_{}'.format(wmp_name[6:-3])

  result.append(nuc_name)
  result.append(wmp_name)
  result.append("[{:e}, {:e}]".format(wmp.start_E, wmp.end_E))
  result.append("{}".format(wmp.data.shape[0]))
  n_wins = int((np.sqrt(wmp.end_E) - np.sqrt(wmp.start_E))/wmp.spacing + 1)
  result.append("{}".format(n_wins))
  result.append("{}".format(wmp.fit_order))
  if len(wmp.pseudo_k0RS) == 0 or (
        len(wmp.pseudo_k0RS) == 1 and wmp.pseudo_k0RS[0] == 0.):
    result.append("Y")
  else:
    result.append(" ")

  nuclides.append(result)

# dump nuclides
output_file = '../Nuclides.md'
format_str = '| {:8} | {:11} | {:28} | {:7} | {:9} | {:8} | {:5} |\n'
table_sep = ['-'*8, '-'*11, '-'*28, '-'*7, '-'*9, '-'*8, '-'*5]
with open(output_file, 'w') as f:
  f.write('# WMP Library Overview\n\n')
  f.write(format_str.format(*headers))
  f.write(format_str.format(*table_sep))
  for nuc in nuclides:
    f.write(format_str.format(*nuc))
