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

- LAMMPS is an open-source software package. System requirements and installation details for LAMMPS can be found on the official LAMMPS website [LAMMPS]([https://www.lammps.org/](https://docs.lammps.org/Install.html).
- The custom C codes used for coarse-grained polymer simulations and analysis were tested on a Linux system (Ubuntu), but should also work on macOS or Windows systems with appropriate compilers.
- Additional data analysis was performed using Python and Jupyter notebooks.
- The C codes were compiled using the GCC compiler with the math library (`-lm`).
