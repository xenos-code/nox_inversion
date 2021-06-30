# nox_inversion

## Files

beta_calc_levs.py
- main program
- contains code to compute NOx emissions inversion from two CMAQ runs

make_intermediate_files.py
- creates model VCD files required by beta_calc_levs
- needs CMAQ CONC files, METCRO2D files

make_emission_files.py
- aggregates NOx emissions data into single files for faster processing
- required by beta_calc_levs

do_inversion.py
- python script to do the inversion and create plots of results
- creates output analysis file used to scale emissions

scale_emis.sh
- scales CMAQ emission files based on analysis file

beta.submit
- script to submit inversion to SLURM on ATMOS

make_intermediate_files.submit
- script to submit intermediate file maker to SLURM on ATMOS


## Steps for inversion
1. make intermediate files (edit file in this dir, submit to SLURM)
2. make emissions files (edit 'emisdir' in make_emissions_files.py in this dir, submit to SLURM)
3. link previous iter analysis file to case dir
4. link previous 2 emis files to case dir
5. edit do_inversion.py in case dir and submit SLURM
6. edit scale_emis.submit in case dir and submit to SLURM 
7. Link emissions files (linking scripts in {emisbasedir}/antbeN_posterior/.
8. Run next iteration

