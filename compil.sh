#!/bin/bash
gfortran precision_mod.f90
gfortran functions_mod.f90
gfortran quad_mod.f90
gfortran quad.f90 functions_mod.f90 quad_mod.f90 precision_mod.f90 -o quad