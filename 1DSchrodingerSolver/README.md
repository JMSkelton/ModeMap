1D Schr&ouml;dinger Equation
============================

A Fortran code for solving the 1D Schr&ouml;dinger equation for a supplied potential.

This code is intended as a post-processing tool for potential-energy surface maps generated using the `ModeMap*` scripts.
In particular, the potentials can be input as a set of polynomial coefficients, which are output by `ModeMap_PolyFit.py`.

The code was written by J. Buckeridge (UCL), and the methods used in the solver are described in [Ref. 1](#Ref1).

Compilation
-----------

The source code needs to be compiled using a Fortran 90 compiler such as GFortran, e.g.:

`gfortran schrodinger-solve1d-anharmonic.f90 -o schrodinger-solve1d-anharmonic.exe`

Usage
-----

The code reads from an input file (`input.dat`), which contains the following information:

1. The number of grid points in the basis set (equal to the number of eigenvalues to find).
2. A scale factor for the effective mass of the phonon - typically 1.0.
3. The solver to use - solution in Fourier space or using the "shooting method".
4. The temperature at which to calculate an effective (renormalised) harmonic frequency.
5. How the potential is input - either from a separate file (`potential.dat`), or as a set of polynomial coefficients.
6. If supplying the potential as a polynomial fit, the number and values of the coefficients and the range of amplitudes over which to sample it.

An example input file with comments explaining the parameters is included with the code (`input.example.dat`).

Once the input file(s) have been set up, the compiled executable can be run with e.g.:

`./schrodinger-solve1d-anharmonic.exe | tee out.dat`

(this will redirect a copy of the output into `out.dat`.)

Convergence
-----------

Some of the input parameters need to be converged, in particular:

* The number of grid points in the basis set; and
* The amplitude range over which the potential/polynomial is evaluated.

For calculating effective harmonic frequencies, convergence with respect to these parameters should be checked at the highest and lowest temperatures to be studied.

Applications
------------

1. *Eigenvalues.*
   After running the program, the calculated eigenvalues (energy levels) are written to `eigenvalue.dat` (in units of eV).
   The energy levels are independent of the supplied temperature, so to obtain them the code only needs to be run once.

2. *Effective (renormalised) harmonic frequencies.*
   This code can calculate effective harmonic frequencies that reproduce the contribution of a mode with an anharmonic potential-energy landscape to the thermodynamic partition function at a given temperature.
   To do this, the solver needs to be called repeatedly with different temperatures set in the input file.
   An example Bash script for running calculations over a range of temperatures is provided for this purpose (`LoopTemperatures.sh`).

Citation
--------

If you use this code in your work, please cite Buckeridge *et al.*, *Phys. Rev. B* **84**, 144120 (**2011**) ([Ref. 1](#Ref1) in the Reference section) and Skelton *et al.*, *Phys. Rev. Lett.* **117**, 075502 (**2016**) ([Ref. 2](#Ref2)).

References
----------

1. <a name="Ref1"></a>J. Buckeridge and S. Fahy, "Mobility in gated GaN<sub>*x*</sub>As<sub>1-*x*</sub> heterostructures as a probe of nitrogen-related electronic states", *Physical Review B* **84**, 144120 (**2016**), DOI: [10.1103/PhysRevB.84.144120](https://doi.org/10.1103/PhysRevB.84.144120)

2. <a name="Ref2"></a>J. M. Skelton, L. A. Burton, S. C. Parker, A. Walsh, C.-E. Kim, A. Soon, J. Buckeridge, A. A. Sokol, C. R. A. Catlow, A. Togo and I. Tanaka, "Anharmonicity in the High-Temperature *Cmcm* Phase of SnSe: Soft Modes and Three-Phonon Interactions", *Physical Review Letters* **117**, 075502 (**2016**), DOI: [10.1103/PhysRevLett.117.075502](https://doi.org/10.1103/PhysRevLett.117.075502)
