# NO2-PESs
Potential energy surfaces for NO2

**Requirements**

(1) RKHS toolkit : Download from https://github.com/MeuwlyGroup/RKHS

**Compile a particular PES**

A test program file (pes_test.f90) is given and it can be compiled for the 2A' PES as

`gfortran RKHS.f90 NO2-2AP-PES.f90 pes_test.f90`

**Running the executable**

Before running the executable make sure that the asymp.dat, .csv and/or .kernel files for that PES present in the present directory (or change the file path in the fortran program).
