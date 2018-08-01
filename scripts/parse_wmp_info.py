#!/usr/bin/env python3

import os
import glob
import numpy as np

import WMP

wmp_dir = "../WMP_Library" # WMP library PATH
wmp_files = glob.glob(os.path.join(wmp_dir, "*.h5"))
wmp_files.sort()

nuclides = []
headers = [
          'Nuclide',
          'WMP File',
          'Energy Range',
          '# Poles',
          '# Windows',
          'CF Order',
         ]

for i, wmp_library in enumerate(wmp_files):
  result = []
  nuc_wmp = WMP.WindowedMultipole.from_hdf5(wmp_library)

  wmp_name = os.path.basename(wmp_library)
  atomic_number = int(wmp_name[0:3])
  mass_number = int(wmp_name[3:6])
  nuc_name = WMP.ATOMIC_SYMBOL[atomic_number] + str(mass_number)
  if wmp_name[6:-3] is not '':
    nuc_name += '_{}'.format(wmp_name[6:-3])

  result.append(nuc_name)
  result.append(wmp_name)
  result.append("[{:e}, {:e}]".format(nuc_wmp.start_E, nuc_wmp.end_E))
  result.append("{}".format(nuc_wmp.data.shape[0]))
  n_wins = int((np.sqrt(nuc_wmp.end_E) - np.sqrt(nuc_wmp.start_E))/nuc_wmp.spacing + 1)
  result.append("{}".format(n_wins))
  result.append("{}".format(nuc_wmp.fit_order))

  nuclides.append(result)

# dump nuclides
output_file = '../nuclides.md'
format_str = '| {:8} | {:11} | {:28} | {:7} | {:9} | {:8} |\n'
table_sep = ['-'*8, '-'*11, '-'*28, '-'*7, '-'*9, '-'*8]
with open(output_file, 'w') as f:
  f.write('# WMP Library Overview\n\n')
  f.write(format_str.format(*headers))
  f.write(format_str.format(*table_sep))
  for nuc in nuclides:
    f.write(format_str.format(*nuc))
