# Magnetizer #

[![arXiv:1809.10521](http://img.shields.io/badge/arXiv-1809.03595-B31B1B.svg)](https://arxiv.org/abs/1809.10521)

This code takes the output of the [Galform][GLF] [semi-analytic model of galaxy formation][SAM]
(SAM) and computes radial dependent ISM properties and magnetic field for each
simulated galaxy. The magnetic field is obtained by numerically solving the
galactic dynamo equation throughout history of each galaxy.

For more details about the physical model, we refer the reader to to the original [code paper][CodePaper].

[CodePaper]: https://ui.adsabs.harvard.edu/#abs/2019MNRAS.483.2424R/
[SAM]: https://ui.adsabs.harvard.edu/#abs/2006RPPh...69.3101B/
[GLF]: https://ui.adsabs.harvard.edu/#abs/2000MNRAS.319..168C

## Quick install and run guide ##

Galaxy Magnetizer requires the following libraries to run:

 * [MPI](https://www.open-mpi.org/)
 * [HDF5](https://www.hdfgroup.org/) - compiled with Fortran and parallel support
 * [GSL](https://www.gnu.org/software/gsl/)
 * [FGSL](http://www.lrz.de/services/software/mathematik/gsl/fortran/)

Once they are installed (see [Dependencies](#dependencies) for building
instructions), the code can be compiled by using `make prod` or `make test`, for
a production or test (debugging and backtracing enable) run, respectively.
For building using multiple processors, the command
`make -j <number of processors>` should work.

The code can be run using mpi:
```
mpirun ./Magnetizer.exe <parameters_file>
```
The parameter file must include the path to the HDF5 input file, containing the 
time-evolving galaxy properties.
The input file can be generated using the scripts/prepare_input.py script.
Parameters not specified in the parameters file are set to their default values,
thus the minimal parameters file is
```
&global_pars
  input_file_name = "sam_input.hdf5"
/
```
An example parameters file can be found in the
`example/example_global_parameters.in`.
In the same directory, there is an example input file,
`example/example_SAM_input.hdf5`.

The Magnetizer can also be run in the _single galaxy mode_ (which is
particularly useful for debugging) by simply specifying the galaxy number i.e.
```
./Magnetizer.exe <parameters_file> <igal> [-f]
```
note that a full output file will be generated, but containing only this galaxy.
The option `-f` allows one to force re-run a particular galaxy.

Magnetizer comes with a range of Python modules and scripts which can be used for

 * preparing the HDF5 input files
 * diagnostic of finished runs
 * visualization and analysis tasks

The python code will depend on the following

 * [numpy](http://www.numpy.org/)
 * [matplolib](http://matplotlib.org/)
 * [h5py](http://www.h5py.org/)


### Preparing an input files ###

Input files can be generated from a Galform (and later other SAMs/sims) run
using the prepare_input.py. For information about this command's
usage/arguments, please check
```
./python/prepare_input.py --help
```
An example Galform output can be found at `scripts/test_SAM_output/galaxies.hdf5`.


### MagnetizerRun objects ###

To facilitate any analysis tasks, one can use the `MagnetizerRun` python object to interact with the input and output. These objects will only load information of galaxies which were completed (i.e. galaxies which did not run either because of error or because Magnetizer was prematurely interrupted are automatically ignored).

```python
run = magnetizer.MagnetizerRun(input_path='example/example_input.hdf5', 
                               output_path='example/example_output.hdf5')
```
It is also possible to supply lists of filenames to input_path and output_path (for example in the case where one is working with various subvolumes of a large simulation, each of them stored in a separate file). The MagnetizerRun object will concatenate the data of the different files in the lists whenever this is needed.

 MagnetizerRun objects contain several useful attributes, as 
```python
zs = run.redshifts
times = run.times
number_of_galaxies = run.ngals
number_of_grid_points = run.ngrid
```

Using this object, for example, one can load the azimuthal component of the field, $B_\phi$, at $z=0$, for all the galaxies using:
```python
Bp_all_galaxies = run.get('Bp', 0.0)
```
where `Bp_all_galaxies` will be a (number_of_galaxies)x(number_of_grid_points) array.

One can also load all the data (i.e. from all redshifts) for a specific galaxy using
```python
Bp_gal_answer = run.get_galaxy('Bp', 42)
```
in this case `Bp_gal_answer` will be a (number_of_grid_points)x(number_of_redshifts) array.


### Contact ###

Please use Github's [tools][issues]
or email [Luiz Felippe S. Rodrigues](mailto:luiz.rodrigues@ncl.ac.uk) if you find any problem.

[issues]: https://github.com/luizfelippesr/magnetizer/issues
