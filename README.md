# stzpp
This repository contains code to reproduce the results of publications [1] and [2].

In addition, the code comprises an extensive simulation packag for high-performance, parallel, MPI-Based simulation of plastic deformation in three-dimensional amorphous glasses via the shear transformation zone theory of amorphous plasticity. While originally used to study deformation of metallic glasses, the software has applications in the simulation of many other forms of disordered materials.

The folder ``mg_3d_template`` also contains a highlight efficient, self-contained, MPI-parallelized implementation of a three-dimensional geometric multigrid method written with ``C++`` templates. This will be released as a standalone package in the coming future.

# Installation
The main simulation packages requires ``mpi`` and ``openmp``. Both can be installed via standard package managers.

# Usage
The main simulation code can be compiled by running ``make``.

The executable ``shear_sim`` can be used to run a simulation of a glassy material undergoing shear deformation between two rigid parallel plates, as well as a simulation of a simple numerical model of friction welding.

The executable ``trans_sim`` can be used to run a simulation of a glassy material subject to an external deformation specified as an abstract coordinate transformation. ``trans_sim.hh`` contains implementations of pure shear, simple shear, and homogeneous deformation. Other transformations can easily be implemented by specifying the corresponding matrix ``T(t)`` and its derivatives in ``trans_sim.hh``.

# Referencing
If you found this repository useful in your research, please consider citing

[1] N. M. Boffi, Chris H. Rycroft, “Parallel three-dimensional simulations of quasi-static elastoplastic solids,” Computer Physics Communications 257, 107254 (2020).

[2] N. M. Boffi, Chris H. Rycroft, “Coordinate transformation methodology for simulating quasi-static elastoplastic solids,” Physical Review E 101, 053304 (2020).

```
@article{boffi_parallel_2020,
	title = {Parallel three-dimensional simulations of quasi-static elastoplastic solids},
	volume = {257},
	issn = {0010-4655},
	doi = {10.1016/j.cpc.2020.107254},
	journal = {Computer Physics Communications},
	author = {Boffi, Nicholas M. and Rycroft, Chris H.},
	month = dec,
	year = {2020},
	pages = {107254},
}
```

```
@article{boffi_coordinate_2020,
	title = {Coordinate transformation methodology for simulating quasistatic elastoplastic solids},
	volume = {101},
	issn = {2470-0045, 2470-0053},
	doi = {10.1103/PhysRevE.101.053304},
	number = {5},
	journal = {Physical Review E},
	author = {Boffi, Nicholas M. and Rycroft, Chris H.},
	month = may,
	year = {2020},
	pages = {053304},
}
```
