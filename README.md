# Windowed Multipole Library

This repository contains neutron cross section data for isotopes in
[windowed multipole (WMP)] format which requires less memory than pointwise
cross sections and allows on-the-fly Doppler broadening to arbitrary temperature.

In WMP format, the entire energy range is chopped into equal-in-momentum domains
, or windows. For each window, poles in certain range are evaluated exactly and
the rest is curve-fitted by low-order polynomials.
The 0K cross section can be expressed as:

![wmp 0k](https://latex.codecogs.com/gif.latex?\sigma(E,&space;T=0\text{K})&space;=&space;\frac{1}{E}&space;\Re\left[\sum_{j=w_1}^{w_2}&space;\frac{ir_j}{\sqrt{E}-p_j}\right]&space;&plus;&space;\sum_{n=0}^{N}&space;c_n&space;E^{\frac{n}{2}-1})

where ![wmp 0k](https://latex.codecogs.com/gif.latex?p_j) and
![wmp 0k](https://latex.codecogs.com/gif.latex?r_j) are poles and residues in the complex
plane and ![wmp 0k](https://latex.codecogs.com/gif.latex?c_n) is the polynomial coefficients.

Both the poles term and polynomial term can be analytically Doppler broadened to
any temperature.

![wmp Tk](https://latex.codecogs.com/gif.latex?\sigma(E,&space;T)&space;\approx&space;\frac{\sqrt{\pi}}{2&space;E&space;\sqrt{\xi}}&space;\sum_j&space;\text{Re}&space;\left[r_j&space;W_i(z)\right]&space;&plus;&space;\sum_{n=0}^{N}c_n\mathcal{D}_n)

![wmp Tk cont](https://latex.codecogs.com/gif.latex?W_i(z)&space;=&space;\frac{i}{\pi}&space;\int_{-\infty}^\infty&space;dt&space;\frac{e^{-t^2}}{z&space;-&space;t},)
![wmp Tk cont](https://latex.codecogs.com/gif.latex?z&space;=&space;\frac{\sqrt{E}&space;-&space;p_j}{2&space;\sqrt{\xi}},)
![wmp Tk cont](https://latex.codecogs.com/gif.latex?\xi&space;=&space;\frac{k_B&space;T}{4&space;A})

![wmp Tk cont](https://latex.codecogs.com/gif.latex?\mathcal{D}_0&space;=&space;\frac{1}{E}erf(\sqrt{\alpha&space;E}),)
![wmp Tk cont](https://latex.codecogs.com/gif.latex?\mathcal{D}_1&space;=&space;\frac{1}{\sqrt{E}},)
![wmp Tk cont](https://latex.codecogs.com/gif.latex?\mathcal{D}_2=&space;\left[\frac{1}{\sqrt{2\alpha}}&space;&plus;&space;E&space;\right]\mathcal{D}_0&space;&plus;&space;\frac{e^{-\alpha&space;E}}{\sqrt{\alpha&space;\pi&space;E}})

![wmp Tk cont](https://latex.codecogs.com/gif.latex?\mathcal{D}_{n&plus;2}=\left[\frac{2n&plus;1}{2\alpha}&space;&plus;&space;E&space;\right]\mathcal{D}_n&space;-&space;\frac{n(n-1)}{4\alpha^2}\mathcal{D}_{n-2},&space;n>0)

The poles of the nuclides in this library have been generated using Vector
Fitting technique [1].
Then the [OpenW] code is used to process the multipole cross sections into
windowed multipole library.

The current library includes 423 nuclides processed from [ENDF/B VII.1] library,
with a target maximum relative error for all cross sections of 0.1%.
An overview of the library can be found in the table in
'[nuclides.md](./nuclides.md)'.

More details about the windowed multipole method or multipole representation can
be found in [2-4].

- [1] J. Liang, X. Peng, S. Liu, C. Josey, B. Forget and K. Smith. "Processing
  Of A Comprehensive Windowed Multipole Library via Vector Fitting." PHYSOR
  2018: Reactor Physics paving the way towards more efficient systems. Cancun,
  Mexico, April 22-26, 2018.
- [2] C. Josey, P. Ducru, B. Forget, and K. Smith. "Windowed multipole for cross
  section Doppler broadening." Journal of Computational Physics, volume 307,
  pp. 715–727 (2016).
- [3] B. Forget, S. Xu, and K. Smith. "Direct Doppler broadening in Monte Carlo
  simulations using the multipole representation." Annals of Nuclear Energy,
  volume 64, pp. 78–85 (2014).
- [4] R. N. Hwang. "A rigorous pole representation of multilevel cross sections
  and its practical applications." Nuclear Science and Engineering, volume
  96(3), pp. 192–209 (1987).

## Library format specifications

Windowed multipole data is stored in HDF5 files, containing the energy
boundaries, window structure, poles/residues, and curve fit coefficients, etc.

Detailed specifications can be found in
'[Windowed Multipole Library Format](./wmp_format.md)'.

## Download

[Git LFS] is used to store the binary HDF5 files.
To download the real WMP library, you need to firstly install [Git LFS] and then
clone this repository.

``` bash
$ git clone https://github.com/mit-crpg/WMP_Library.git
```

## Usage

With the windowed multipole library, it is convenient and efficient to compute
cross sections at arbitrary temperature (in 0 K-3000 K range) and energy
(in the resolved resonance range) for 4 reactions: total, elastic scattering,
absorption, and fission.

An excellent reference is [OpenMC], which implements WMP in both the transport
solver and the [OpenMC Python API]. You can also check the scripts `scripts/WMP.py`
for an Python implementation, which illustrates how to read, evaluate and export
an windowed multipole format library.
For example, the following script demonstrates how to utilize WMP library using
the nuclear data interface [WindowedMultipole].

``` python
# load a library
import WMP
u238_multipole = WMP.WindowedMultipole.from_hdf5('092238.h5')

# evaluate cross sections at a given energy and temperature
scatt_xs, absorption_xs, fission_xs = u238_multipole(E=1.0, T=300.)

# comparison with ACE library (HDF5 format used in OpenMC)
import openmc.data
u238_ace = openmc.data.IncidentNeutron.from_hdf5('U238.h5')
energy = np.logspace(np.log10(u238_multipole.E_min), np.log10(u238_multipole.E_max), 1E4)
scatt_xs_wmp = u238_multipole(energy, T=293.75)[0]
scatt_xs_ace = u238_ace[2].xs['294K'](energy)
# then you can plot the cross sections with energies
```

## Reporting

 - Submit GitHub issues: https://github.com/mit-crpg/WMP_Library/issues
 - Contact MIT [Computational Reactor Physics Group (CRPG)]

[windowed multipole (WMP)]: http://openmc.readthedocs.io/en/latest/methods/cross_sections.html#windowed-multipole-representation
[OpenMC]: https://github.com/mit-crpg/openmc
[OpenW]: https://github.com/mit-crpg/WHOPPER
[ENDF/B VII.1]: http://www.nndc.bnl.gov/endf/b7.1/
[Windowed Multipole Library Format]: http://openmc.readthedocs.io/en/latest/io_formats/data_wmp.html#io-data-wmp
[Git LFS]: https://git-lfs.github.com
[OpenMC Python API]: http://openmc.readthedocs.io/en/latest/pythonapi/index.html
[WindowedMultipole]: https://github.com/mit-crpg/WMP_Library/blob/master/scripts/WMP.py
[Computational Reactor Physics Group (CRPG)]: http://crpg.mit.edu/
