[![DOI](https://zenodo.org/badge/364278274.svg)](https://zenodo.org/badge/latestdoi/364278274)

# SPECFEM2D-DG

SPECFEM2D-DG allows users to perform 2D simulations of wave propagation in the mechanically-coupled solid-atmosphere system.
The solid module resolves the equations of visco-elastodynamics; it is a direct import from the well-known [SPECFEM2D](https://github.com/geodynamics/specfem2d) software.
The atmospheric modules resolve the Navier-Stokes equations, either the classical non-linear ones or linearised ones.

Users interested in modelling wave propagation in linear acoustic media (incompressible windless non-viscous air or water) or poroelastic media are referred to [SPECFEM2D](https://github.com/geodynamics/specfem2d) and [SPECFEM3D](https://github.com/geodynamics/specfem3d).

The main "historical" developers of SPECFEM are Dimitri Komatitsch and Jeroen Tromp, though there are currently many more.
The main "historical" developers of the DG extension are Quentin Brissaud and Léo Martire.


## Configuration and Installation

Refer to the `README_SPECFEM2D_DG.md` file.


## Development

Development is hosted on GitHub in [project SAMoSA's specfem2d-dg repository](https://github.com/samosa-project/specfem2d-dg).


## Citation

- Full Navier-Stokes module:
_Q. Brissaud, R. Martin, R. F. Garcia, and D. Komatitsch, “[Hybrid Galerkin numerical modelling of elastodynamics and compressible Navier-Stokes couplings: Applications to seismo-gravito acoustic waves](https://academic.oup.com/gji/article/210/2/1047/3798201)”, Geophys. J. Int., vol. 210, no. 2, pp. 1047–1069, 2017, doi: 10.1093/gji/ggx185._
- Linear Navier-Stokes module:
_L. Martire, R. Martin, Q. Brissaud, and R. F. Garcia, “[SPECFEM2D-DG, an Open Source Software Modeling Mechanical Waves in Coupled Solid-Fluid Systems: the Linearised Navier-Stokes Approach](https://academic.oup.com/gji/advance-article/doi/10.1093/gji/ggab308/6342174?guestAccessKey=3c876f08-2039-4adc-8c7d-9173b208d1c8)”, Geophys. J. Int., 2021, doi: 10.1093/gji/ggab308._
```
@article{Brissaud2017_SPECFEM2D_DG_FNS,
author = {Brissaud, Quentin and Martin, Roland and Garcia, Rapha{\"{e}}l F. and Komatitsch, Dimitri},
doi = {10.1093/gji/ggx185},
journal = {Geophysical Journal International},
number = {2},
pages = {1047--1069},
title = {{Hybrid Galerkin numerical modelling of elastodynamics and compressible Navier-Stokes couplings: Applications to seismo-gravito acoustic waves}},
volume = {210},
year = {2017}
}
@article{Martire2021_SPECFEM2D_DG_LNS,
author = {Martire, L{\'{e}}o and Martin, Roland and Brissaud, Quentin and Garcia, Rapha{\"{e}}l F.},
doi = {10.1093/gji/ggab308},
journal = {Geophysical Journal International},
title = {{SPECFEM2D-DG, an Open Source Software Modeling Mechanical Waves in Coupled Solid-Fluid Systems: the Linearised Navier-Stokes Approach}},
year = {2021}
}
```
