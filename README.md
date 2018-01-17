# Magnetizer #

This code takes the output of a [semi-analytic model of galaxy formation][SAM]
(SAM) or cosmological hydro-simulation _(planned)_ and computes radial dependent
ISM properties and magnetic field for each simulated galaxy. The magnetic field
is obtained by numerically solving the galactic dynamo equation throughout
history of each galaxy.

[SAM]: https://ui.adsabs.harvard.edu/#abs/2006RPPh...69.3101B/

## Supported SAMs and simulations ##

For the present version, only the [Galform][GLF] model is supported.

We have plans to include suport to the [Galacticus][GLC]. We also have plans to
allow using data from the [Millenium database][MIL], which includes the output
of several published versions of both Galform and also various versions of the
[L-Galaxies][LGA] SAM developed at the MPA. Similarly, there is planned support
to retrieving data from the [Eagle database][EAG], which contains data extracted
from the [EAGLE simulation][EAGD].

[GLF]: https://ui.adsabs.harvard.edu/#abs/2000MNRAS.319..168C
[GLC]: https://sites.google.com/site/galacticusmodel/
[MIL]: http://wwwmpa.mpa-garching.mpg.de/millennium/#DATABASE_ACCESS
[LGA]: http://galformod.mpa-garching.mpg.de/public/LGalaxies/
[EAG]: http://icc.dur.ac.uk/Eagle/database.php
[EAGD]: http://icc.dur.ac.uk/Eagle/

## Quick install and run guide ##

Galaxy Magnetizer requires the following libraries to run:

 * [MPI](https://www.open-mpi.org/)
 * [HDF5](https://www.hdfgroup.org/) - compiled with Fortran and parallel support
 * [GSL](https://www.gnu.org/software/gsl/)
 * [FGSL](http://www.lrz.de/services/software/mathematik/gsl/fortran/)

Once they are installed (see [Dependencies](#dependencies) for building
instructions), the code can be compiled by simply typing `make`
(or `make -j <number of processors>` if you would like to save time using
multiple processors) and run using:
```
mpirun ./Magnetizer.exe <parameters_file>
```
The parameter file (more details in the dedicated section below) must include
the path to the HDF5 input file, containing the time-evolving galaxy properties.
The input file can be generated using the scripts/prepare_input.py script.
Parameters not specified in the parameters file are set to their default values,
thus the minimal parameters file is
```
&global_pars
  input_file_name = "sam_input.hdf5"
/
```
An example parameters file can be found in the `example/example_global_parameters.in`
(actually, if one tries to run the code without specifying a parameter file,
this is the fall-back case). In the same directory, there is an example input
file, `example/example_SAM_input.hdf5`.

The magnetizer can also be run in the _single galaxy mode_ (which is
particularly useful for debugging) by simply specifying the galaxy number i.e.
```
./magnetize_galform.exe <parameters_file> <igal>
```
note that the output file will contain only this galaxy.

Magnetizer comes with a range of Python modules and scripts which can be used for

 * preparing the HDF5 input files
 * diagnostic of finished runs
 * visualization and analysis tasks

The python code will depend on the following

 * [numpy](http://www.numpy.org/)
 * [matplolib](http://matplotlib.org/)
 * [h5py](http://www.h5py.org/)


### Depencencies ###

* Building details

### Preparing an input files ###

Input files can be generated from a Galform (and later other SAMs/sims) run
using the prepare_input.py. For information about this command's
usage/arguments, please check
```
./scripts/prepare_input.py --help
```
An example Galform output can be found at `scripts/test_SAM_output/galaxies.hdf5`.

#### Synthetic SAM output files ####

If one is interested in generating a synthetic example, usually for the purposes
of *testing*, it is possible to do it using the command
```
./scripts/generate_galform_output_from_txt.py <your_galaxy_properties_table.txt
```
Again, run it with `--help` to get usage instructions.

Note that this command assumes each galaxy experiences no mergers, no accretion
and no changes in structure or rotation. The only changes it accounts for are in
the mass of gas: the present day SFR is assumed to be constant and the stellar
masses and gas masses change in time according to it.

An example galaxy properties table containing data from M31 and Milky Way can be
found in `scripts/example_galaxies.txt`.


### Parameter files ###

* Description of how to use the parameter files

### Input files ###


### Scripts ###

* List of the accompanying scripts and python modules
* Code review
* Other guidelines

### Contact ###

Please use Bitbucket's [tools][issues]
or email [Luiz Felippe S. Rodrigues](mailto:luiz.rodrigues@ncl.ac.uk) if you find any problem.

[issues]: https://bitbucket.org/luizfelippe/magnetic-fields-and-galaxy-formation-models/issues
