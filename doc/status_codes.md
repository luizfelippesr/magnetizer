
#       Exit status codes       #

For each galaxy processed and for each snapshot, Magnetizer will
store a 1-character _status code_, which can be
tell whether computation was successful and, if not, what happened.

Here is a list of the possible status codes.


- `M` - Everything went well, both computation of the galaxy profiles and the
        solution of the dynamo equations

- `m` - Everything went well, both computation of the galaxy profiles and the
        solution of the dynamo equations. However, the timestep had to be
        reduced.

- `t` - The computation of galaxy profiles went well. However, the paramter
        p_no_magnetic_fields_test_run was set to T.

- `e` - Disk gas mass or disk size were below thresholds set by
        (possibly, the galaxy became an elliptical from one snapshot to the
        other).

- `g` - Gas density became lower than threshold.

- `h` - Scale-height at the half-mass radius became larger
        p_height_tolerance times the half-mass radius. (I.e. the disc is not
        thin!)

- `P` - Invalid solution to the hydrotatic equilibrium equation detected for a
        one specific radius. A patch was applied: the neighbouring values for
        h and rho were used instead.
        NB This is usually harmless, if the patch is applied multiple times the code will fail with status H instead of this.

- `H` - Negative scale-height detected. Most likely, the code failed to solve
        the vertical pressure equilibrium and the galaxy aborted.

- `s` - Timestep became too short and the galaxy was aborted.

- `i` - After reducing the timestep multiple times, the run was aborted.

- `p` - Error while reading input parameters.

- `-` - Redshift without have any data (it was absent in the input files).






