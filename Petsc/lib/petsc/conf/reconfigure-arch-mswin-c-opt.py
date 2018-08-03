#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--prefix=/home/alireza/PetscInstall',
    '--with-blaslapack-dir=/cygdrive/E/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64',
    '--with-cc=win32fe icl',
    '--with-debugging=0',
    '--with-fc=win32fe ifort',
    '--with-mpi-include=/cygdrive/E/Program Files (x86)/IntelSWTools/mpi/2018.1.156/intel64/include',
    '--with-mpi-lib=[/cygdrive/E/Program Files (x86)/IntelSWTools/mpi/2018.1.156/intel64/lib/release/impi.lib]',
    '--with-mpi-mpiexec=/cygdrive/E/Program Files (x86)/IntelSWTools/mpi/2018.1.156/intel64/bin/mpiexec.exe',
    '-CFLAGS=-O2 -MT -wd4996 -Qopenmp',
    '-CXXFLAGS=-O2 -MT -wd4996 -Qopenmp',
    '-FFLAGS=-MT -O2 -Qopenmp',
    'PETSC_ARCH=arch-mswin-c-opt',
    '–with-openmp',
  ]
  configure.petsc_configure(configure_options)
