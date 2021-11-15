
#Following commands can be used to run observables_singles


./Observables_single.exe "RM" files_luiz/Lacey14_fiducial_1.in 0 40 0.7 0.002 0  0.1 0.1 /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/
WHAT IS THE RANGE OF X AND Y AND IN GENERAL OTHER PARAMETERS?

./Observables_single.exe "RM" files_luiz/Lacey14_fiducial_1.in 9 40 0.7 0.002 0  /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/

./Observables_single.exe "I" /home/darkmatter/work/charles/files_luiz/Lacey14_fiducial_1.in 1 12 0.7 0.002 0 0 

./Observables_single.exe "Image" files_luiz/Lacey14_fiducial_1.in 9 40 0.7 0.002 0  /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/

./Observables_single.exe "RM_study" files_luiz/Lacey14_fiducial_1.in 9 40 0.7 0.002 0  /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/

./Observables_single.exe "RM_study" files_luiz/Lacey14_fiducial_1.in 9 1 0.7 0.002 0  /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/
 

PRODUCED AN OUTPUT. BUT DIDN'T PRITNT ANY

./Observables_single.exe "FRB" files_luiz/Lacey14_fiducial_1.in 9 40 0.7 0.002 0  /home/carlos/work/magnetic_fields/magnetizer-master/files_luiz/



# case with lot of error

./Observables_single.exe "I" example/example_global_parameters.in 10 12 0.7 0.002 0 0
./Observables.exe 'example/example_global_parameters.in' 'synchrotron' 0.0
