Data-driven polymer modeling reveals how scale-dependent active fluctuations shape chromatin organization

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [File Description](#file-description)
- [How to Run](#how-to-run)

## Overview

This repository contains the simulation and analysis codes used in the study:

"Data-driven Polymer Modeling Reveals How Scale-Dependent Active Fluctuations Shape Chromatin Organization"

We use the open-source LAMMPS package to simulate chromatin at nucleosome resolution. Custom C codes are used for simulating coarse-grained chromatin polymers. All analyses were performed using C and Python scripts.

## System Requirements

- LAMMPS is an open-source software package. System requirements and installation details for LAMMPS can be found on the official LAMMPS website [LAMMPS]([https://www.lammps.org/](https://docs.lammps.org/Install.html)).
- The custom C codes used for coarse-grained polymer simulations and analysis were tested on a Linux system (Ubuntu), but should also work on macOS or Windows systems with appropriate compilers.
- Additional data analysis was performed using Python and Jupyter notebooks.
- The C codes were compiled using the GCC compiler with the math library (`-lm`).

## File Description

- `run_microc_simulation.in` - Script file for LAMMPS
- `cg_analysis.c` - Code to measure coarse-graining properties
- `CG_activity_simulation.c` - Code for simulating coarse-grained chromatin polymers with local activities (LA)
- `CG_activity_with_hic_simulation.c` - Code for simulating coarse-grained chromatin polymers with local activities (LA) and Hi-C-based interactions
- `fit.py` - Code for finding the best linear combination of distributions that fits the target distribution


## How to Run

Refer to the official [LAMMPS installation documentation](https://docs.lammps.org/Install.html) for installation instructions.
Run the nucleosome-resolution chromatin simulation using:
```bash
./lmp_serial < run_microc_simulation.in
```
Compile and run `cg_analysis.c` using:
```bash
gcc cg_analysis.c -lm
./a.out
```
Compile and run `CG_activity_simulation.c` using:
```bash
gcc CG_activity_simulation.c -lm
./a.out
```
Compile and run `CG_activity_with_hic_simulation.c` using:
```bash
gcc CG_activity_with_hic_simulation.c -lm
./a.out
```
Run the Python analysis script using:
```bash
python fit.py
```
---

**## Expected Output

The simulation generates single position trajectory file. The analysis codes compute coarse-grained polymer properties, while the Python fitting script estimates the parameters that best reproduce the target distributions.
---**

## Typical Runtime

Typical installation time for LAMMPS is around 30 minutes. A single nucleosome-resolution simulation may take several hours depending on the system size and hardware configuration. The coarse-grained simulations also typically require several hours depending on the system size and simulation parameters.
