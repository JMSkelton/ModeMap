ModeMap
=======

A set of tools for mapping and analysing potential-energy surfaces along phonon modes.

<hr>

**Note (26/04/2017):**

*A bug was recently identified in the `ModeMap.py` script whereby normal-mode coordinates were scaled internally by 1/sqrt(N<sub>a</sub>), leading to U(Q) curves that were inconsistent with the harmonic mode frequencies.*
*This issue has now been fixed.*

*If you have output generated with earlier versions of the code, you can use the `--rescale_na` command-line option to the `ModeMap_PostProcess.py` script to correct for the bug during post processing.*
*If you need to check whether you were using the affected version, you can type `ModeMap_PostProcess.py -h` and see whether the `--rescale_na` option appears in the list of command-line arguments.*

*You can also use the new `--harmonic_fit` option of `ModeMap_PolyFit.py` to fit a potential to a harmonic function U(Q) = 1/2 &omega;<sup>2</sup> Q<sup>2</sup> and extract the fitted frequency.*

*We apologise for any inconvenience this may cause, and we are very grateful to Janine George (RWTH Aachen) for pointing out this issue and suggesting a fix.*

<hr>

Prerequisites
-------------

These tools are based around the [Phonopy code](https://atztogo.github.io/phonopy/).
To use them requires a harmonic phonon calculation with Phonopy to have been performed on the system of interest.
It is assumed that users are comfortable with these calculations, and with the basics of the underlying theory.

The `ModeMap*` scripts are written in Python.
They were developed and tested on Python 3.x, but should be compatible with Python >= 2.7.
`ModeMap.py` uses routines from the `Phonopy` Python package to generate structures modulated along phonon modes.
`Phonopy` therefore needs to be installed and "importable" (e.g. the library directories added to `PYTHONPATH`).
`ModeMap_PostProcess.py` and `ModeMap_PolyFit.py` require the `NumPy`, `SciPy` and `Matplotlib` packages.

Users wishing to use the 1D Schr&ouml;dinger solver provided for post processing will need to compile it using a Fortran 90 compiler such as GFortran.

Usage
-----

A typical usage will consist of three or four steps, implemented in a series of short programs:

1. Prepare a sequence of structures with displacements along one (1D map) or two (2D map) phonon modes over a range of amplitudes "frozen in" (`ModeMap.py`).

2. Run single-point energy calculations on these structures, and extract the total energies.
   A basic script for the Vienna *ab initio* Simulation Package (VASP) code is provided (`ExtractTotalEnergies.py`).

3. Post-process the calculation results to generate the potential-energy surface maps (`ModeMap_PostProcess.py`).

4. If desired, fit the calculated potential-energy surfaces to polynomial functions for further analysis (`ModeMap_PolynomialFit.py`).

5. Again, if desired, carry out further post processing by using the fitted potential-surface(s) as input to a 1D Schr&ouml;dinger solver.
   A Fortran code written by J. Buckeridge (UCL) is provided in [1DSchrodingerSolver](./1DSchrodingerSolver).
   This uses the Fourier method to obtain the eigenvalues and eigenvectors in the anharmonic potential, and to determine an effective renormalised (harmonic) frequency for the mode that reproduces its contribution to the thermodynamic partition function at a given temperature.
   The [TISH code](https://github.com/jarvist/Julia-SoftModeTISH-DeformationPotential) by J. M. Frost, which was written to study bandgap-deformation potentials, can also work with polynomial fits.

Limitations
-----------

The expected use case for this code is to study phonon modes along which the potential-energy surface is anharmonic (e.g. modes representing double-well structural instabilities).

By explicitly mapping the potential energy along the mode eigenvector, the energetic minima along imaginary modes can be located.
By sloving a 1D Schr&ouml;dinger equation for the potential, anharmonic eigenvalues (energy levels) can be recovered.
Further analyses can be performed using this as a base, e.g. to study the effect of the soft-mode instabilities on the electronic structure.

Users should bear in mind, however, that this technique makes several assumptions, in particular:

* That the potential-energy surfaces along the phonon modes are independent; and
* That the mode eigenvectors obtained within the harmonic approximation are a reasonable representation of the anharmonic motion.

A common case where these assumptions may not hold is when modes are degenerate.
In such cases, linear combinations of the mode eigenvectors are valid solutions to the harmonic problem.
In the case of degenerate double-well modes, this means that the eigenvectors may not "point" directly at the energetic minima, and the minima instead lie at linear combinations of both eigenvectors.

Within this approximation, the analysis of an *N*-fold degenerate set of modes should technically be approached by mapping the *N*-dimensional potential-energy surface and solving the corresponding *N*-dimensional Schro&ouml;dinger equation.
This is presently beyond the scope of these codes, although they may provide a basis for further development.

Examples
--------

The following examples are provided to illustrate some of the applications of these codes:

* [Cmcm *SnSe*](./Example_Cmcm-SnSe) Reproduces some of the calculations and analysis in [Ref. 1](#Ref1).

Further Information
-------------------

The basic theory and its application to the high-temperature *Cmcm* phase of SnSe are described in detail in [Ref. 1](#Ref1) in the References section.

The mapping codes were used to map the potential-energy surfaces along the imaginary phonon modes in the cubic phase of methylammonium lead iodide ((CH<sub>3</sub>)(NH<sub>3</sub>)PbI<sub>3</sub>, MAPbI<sub>3</sub>) in [Ref. 2](#Ref2).

An extension of the method to study the effect of the phonon instabilities in MAPbI<sub>3</sub> on the electronic structure is described in [Ref. 3](#Ref3).

This method implemented in these codes is very similar in principle to the "Decoupled Anharmonic Mode Approximation" technique developed in [Ref. 4](#Ref4), albeit with the simplification that in this implementation the potential-surface mapping and analysis would normally only be performed for selected modes.

Citation
--------

If you use these codes in your own work, please cite Skelton *et al.*, *Phys. Rev. Lett.* **117**, 075502 (**2016**) ([Ref. 1](#Ref1) in the References section).

References
----------

1. <a name="Ref1"></a>J. M. Skelton, L. A. Burton, S. C. Parker, A. Walsh, C.-E. Kim, A. Soon, J. Buckeridge, A. A. Sokol, C. R. A. Catlow, A. Togo and I. Tanaka, "Anharmonicity in the High-Temperature *Cmcm* Phase of SnSe: Soft Modes and Three-Phonon Interactions", *Physical Review Letters* **117**, 075502 (**2016**), DOI: [10.1103/PhysRevLett.117.075502](https://doi.org/10.1103/PhysRevLett.117.075502)

2. <a name="Ref2"></a>A. N. Beecher, O. E. Semonin, J. M. Skelton, J. M. Frost, M. W. Terban, H. Zhai, A. Alatas, J. S. Owen, A. Walsh and S. J. L. Billinge, "Direct Observation of Dynamic Symmetry Breaking above Room Temperature in Methylammonium Lead Iodide Perovskite", *ACS Energy Letters* **1** (4), 880 (**2016**), DOI: [10.1021/acsenergylett.6b00381](https://doi.org/10.1021/acsenergylett.6b00381)

3. <a name="Ref3"></a>L. D. Whalley, J. M. Skelton, J. M. Frost and A. Walsh, "Phonon anharmonicity, lifetimes, and thermal transport in CH<sub>3</sub>NH<sub>3</sub>PbI<sub>3</sub> from many-body perturbation theory", *Physical Review B* **94**, 220301(R) (**2016**), DOI: [10.1103/PhysRevB.94.220301](https://doi.org/10.1103/PhysRevB.94.220301)

4. <a name="Ref4"></a>D. J. Adams and D. Passerone, "Insight into structural phase transitions from the decoupled anharmonic mode approximation", *Journal of Physics: Condensed Matter* **28** (30), 305401 (**2016**), DOI: [10.1088/0953-8984/28/30/305401](http://dx.doi.org/10.1088/0953-8984/28/30/305401)