1)copy the custom pair style files into the lammps src directory

2)install lammps with cmake:
cd lammps                # change to the LAMMPS distribution directory
mkdir build; cd build    # create and use a build directory
cmake3 ../cmake -D PKG_COLLOID=yes -D PKG_BROWNIAN=yes -D PKG_DIPOLE=yes # configuration reading CMake scripts from ../cmake
cmake3 --build . # compilation (or type "make")
make install

3)run lammps:
lmp -in brownian.lammps -screen none &
or in parallel: mpirun -np 2 lmp -in brownian.lammps -screen none &
