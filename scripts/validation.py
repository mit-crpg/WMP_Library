#!/usr/bin/env python3

import os
import glob
import time
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import openmc.data
import WMP

from optparse import OptionParser

WMP_PATH = "../WMP_Library" # WMP library PATH
ACE_PATH = "../../njoy_293.75K" # ACE library PATH
OUT_PATH = "../WMP_Validation" # OUTPUT PATH
TEMPERATURE = 293.75

# Command line parsing
usage = """usage: %prog [options]"""
parser = OptionParser(usage=usage)
parser.add_option('-w', '--wmp_directory', dest='wmpdir', default=WMP_PATH,
                  help="Directory for windowed multipole library. "
                  "Default: {}".format(WMP_PATH))
parser.add_option('-f', '--wmp_file', dest='wmpfile',
                  help="Specify the wmp file to process. ")
parser.add_option('-a', '--ace_directory', dest='acedir', default=ACE_PATH,
                  help="Directory for ACE library. "
                  "Default: {}".format(ACE_PATH))
parser.add_option('-o', '--out_directory', dest='outdir', default=OUT_PATH,
                  help="Directory for outputs. "
                  "Default: {}".format(OUT_PATH))
parser.add_option('-t', '--temperature', dest='temp', default=TEMPERATURE,
                  type='float', help="Temperature to compare. "
                  "Default: {}".format(TEMPERATURE))
(options, args) = parser.parse_args()

wmp_dir = options.wmpdir
ace_dir = options.acedir
out_dir = options.outdir
temp = options.temp
strTemp = "{}K".format(int(round(temp)))

if options.wmpfile is not None:
  assert os.path.isfile(options.wmpfile), "wmp library {} not found".format(options.wmpfile)
  wmp_files = [options.wmpfile]
else:
  assert os.path.exists(wmp_dir), "wmp library dir {} not found".format(wmp_dir)
  wmp_files = glob.glob(os.path.join(wmp_dir, "*.h5"))

assert os.path.exists(ace_dir), "ace library dir {} not found".format(ace_dir)

if not os.path.exists(out_dir):
  os.makedirs(out_dir)

print("Start validating {} nuclides - {}".format(len(wmp_files), time.ctime()))
for i, wmp_library in enumerate(wmp_files):
  # load wmp data
  nuc_wmp = WMP.WindowedMultipole.from_hdf5(wmp_library)
  nuc_name = nuc_wmp.name

  print("{:>3}/{:<3} Processing {} {} - {} ".format(
          i+1, len(wmp_files), nuc_name, wmp_library, time.ctime()))

  ace_file = os.path.join(ace_dir, nuc_name+'.h5')
  assert os.path.isfile(ace_file), "ace_file {} not found".format(ace_file)

  logfile_name = '{}_{}K_validation.log'.format(nuc_name, temp)
  logfile = os.path.join(out_dir, logfile_name);
  f = open(logfile, 'w');

  # write info
  f.write("WMP file: {}\n".format(wmp_library))
  f.write("Nuclide: {}\n".format(nuc_name))
  f.write("Energy range: [{}, {}] eV\n".format(nuc_wmp.E_min, nuc_wmp.E_max))
  f.write("Number of windows: {}\n".format(nuc_wmp.windows.shape[0]))
  f.write("Fissionable: {}\n".format(nuc_wmp.fissionable))

  # load ace data
  nuc_ace = openmc.data.IncidentNeutron.from_hdf5(ace_file)
  assert strTemp in nuc_ace.temperatures, "ace file does not contain T={}".format(strTemp)

  f.write("Load ace file: {}\n".format(ace_file))

  # energy grid for comparison
  max_e = nuc_wmp.E_max
  min_e = nuc_wmp.E_min
  N_points = 1E4
  energy = np.logspace(np.log10(min_e), np.log10(max_e), N_points)
  energy[0] = min_e
  energy[-1] = max_e
  f.write("Test energy range: [{}, {}] eV\n".format(energy[0], energy[-1]))
  f.write("Test temperature: {} K\n".format(temp))

  # reactions for comparison
  mts = [1, 2, 27, 18]
  reactions = ['total', 'elastic', 'absorption', 'fission']

  # compute cross sections
  xs_wmp = np.zeros((len(mts), len(energy)))
  xs_ace = np.zeros((len(mts), len(energy)))

  xs_wmp[[1,2,3], :] = nuc_wmp(energy, temp)
  xs_wmp[0, :] = xs_wmp[1, :] + xs_wmp[2, :]
  for i, mt in enumerate(mts):
    if mt in nuc_ace:
      xs_ace[i, :] = nuc_ace[mt].xs[strTemp](energy)

  if np.any(xs_wmp < 0.):
    print("!!! Found negative cross sections in WMP library!")

  # compare
  for i, rxn in enumerate(reactions):
    print("  -- {} cross section".format(rxn))
    rxn_wmp = xs_wmp[i]
    rxn_ace = xs_ace[i]
    error = abs(rxn_wmp - rxn_ace)

    if not rxn_ace.any():
      continue

    # max abs. error
    max_error = max(error)
    max_error_idx = np.argmax(error)
    max_error_energy = energy[max_error_idx]
    max_error_wmpxs = rxn_wmp[max_error_idx]
    max_error_acexs = rxn_ace[max_error_idx]
    f.write("{} - max abs error:\n".format(rxn))
    f.write("  energy: {}\n".format(max_error_energy))
    f.write("  WMP xs: {}\n".format(max_error_wmpxs))
    f.write("  ACE xs: {}\n".format(max_error_acexs))
    f.write("  error : {}\n".format(max_error))

    # max rel. error
    relerr = abs(rxn_wmp/rxn_ace - 1)
    relerr[rxn_ace == 0] = 0
    relerr2 = np.array(relerr)
    relerr2[error <= 1E-5] = 0
    max_error = max(relerr2)
    max_error_idx = np.argmax(relerr2)
    max_error_energy = energy[max_error_idx]
    max_error_wmpxs = rxn_wmp[max_error_idx]
    max_error_acexs = rxn_ace[max_error_idx]
    f.write("{} - max rel error:\n".format(rxn))
    f.write("  energy: {}\n".format(max_error_energy))
    f.write("  WMP xs: {}\n".format(max_error_wmpxs))
    f.write("  ACE xs: {}\n".format(max_error_acexs))
    f.write("  error : {:.2f}%\n".format(max_error*100))

    # plot
    plt.clf()
    fig, ax1 = plt.subplots()
    lns1 = ax1.loglog(energy, rxn_wmp, 'g', label="WMP xs")
    lns2 = ax1.loglog(energy, rxn_ace, 'b', label="ACE xs")
    ax2 = ax1.twinx()
    lns3 = ax2.loglog(energy, relerr, 'r', label="rel. err.", alpha=0.5)
    lns = lns1 + lns2 + lns3
    labels = [l.get_label() for l in lns]
    ax1.legend(lns, labels, loc='best')
    ax1.set_xlabel('energy (eV)')
    ax1.set_ylabel('cross section (b)', color='b')
    ax1.tick_params('y', colors='b')
    ax2.set_ylabel('relative error', color='r')
    ax2.tick_params('y', colors='r')

    plt.title("{} {} xs {}K".format(nuc_name, rxn, temp))
    fig.tight_layout()
    figfile = os.path.join(out_dir, "{}_validation_{}K_{}.png".format(nuc_name, temp, rxn))
    plt.savefig(figfile, dpi=600)
    plt.close()

  f.close()

print("Done! - {}".format(time.ctime()))
