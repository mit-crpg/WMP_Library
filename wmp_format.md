# Windowed Multipole Library Format

**/version** (*char[]*)

  The format version of the file.  The current version is "v1.0"

**/nuclide/**

- **broaden_poly** (*int[]*)

    If 1, Doppler broaden curve fit for window with corresponding index.
    If 0, do not.

- **curvefit** (*double[][][]*)

    Curve fit coefficients. Indexed by (reaction type, coefficient index,
    window index).

- **data** (*complex[][]*)

    Complex poles and residues. Each pole has a corresponding set of
    residues. For example, the i-th pole and corresponding residues
    are stored as

   ![data](https://latex.codecogs.com/gif.latex?\text{data}[:,i]&space;=&space;[\text{pole}_i,&space;~\text{residue}_{i1},&space;~\text{residue}_{i2},&space;~\text{residue}_{i3}])

    The residues are in the order: scattering, absorption, fission. Complex
    numbers are stored by forming a type with "r" and "i" identifiers,
    similar to how [h5py] does it.

- **E_max** (*double*)

    Highest energy the windowed multipole part of the library is valid for.

- **E_min** (*double*)

    Lowest energy the windowed multipole part of the library is valid for.

- **spacing** (*double*)

   ![spacing](https://latex.codecogs.com/gif.latex?\frac{\sqrt{E_{max}}-&space;\sqrt{E_{min}}}{n_w})

    Where ![data](https://latex.codecogs.com/gif.latex?E_{max}) is the
    maximum energy the windows go up to.
    ![data](https://latex.codecogs.com/gif.latex?E_{min}) is the minimum
    energy and equivalent to ``E_min``, and ![data](https://latex.codecogs.com/gif.latex?n_w)
    is the number of windows, given by ``windows``.

- **sqrtAWR** (*double*)

    Square root of the atomic weight ratio.

- **windows** (*int[][]*)

    The poles to start from and end at for each window. windows[i, 0] and
    windows[i, 1] are, respectively, the indexes (1-based) of the first and last
    pole in window i.

[h5py]: http://docs.h5py.org/en/latest/
