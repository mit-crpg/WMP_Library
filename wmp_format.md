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

- **end_E** (*double*)

    Highest energy the windowed multipole part of the library is valid for.

- **spacing** (*double*)

   ![spacing](https://latex.codecogs.com/gif.latex?\frac{\sqrt{E_{max}}-&space;\sqrt{E_{min}}}{n_w})

    Where ![data](https://latex.codecogs.com/gif.latex?E_{max}) is the
    maximum energy the windows go up to.
    ![data](https://latex.codecogs.com/gif.latex?E_{min}) is the minimum
    energy and equivalent to ``start_E``, and ![data](https://latex.codecogs.com/gif.latex?n_w)
    is the number of windows, given by ``windows``.

- **sqrtAWR** (*double*)

    Square root of the atomic weight ratio.

- **start_E** (*double*)

    Lowest energy the windowed multipole part of the library is valid for.

- **windows** (*int[][]*)

    The pole to start from for each window. windows[i, 0] - 1 and windows[i, 1]
    -1  are the indexes of the first and last pole in window i respectively.

[h5py]: http://docs.h5py.org/en/latest/
