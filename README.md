# laser_heat
laser ablation_modeling
This is the code for paper "Effect of silicon target porosity on laser ablation threshold: molecular dynamics simulation"
jupyter_create_porous_Si script can be used for generation of the position of silicon atoms in porous silicon. The pore distribution as well as pores shape could be adjusted. Script yields *.csv files (like silicon2_5nmpores0.3.csv in the same folder) with coordinates of all atoms in wafer. Script requires Si_1010100.csv file at the same folder as initial atom positioning in crystalline silicon
heat_source.ipynb is a script for calculation of heat distribution using one-photon and two-photon absorbance coefficients (which lies into silicon_1photon_extinction.txt  and silicon_2photon_absorbance.txt). Missing values are filled using interpolation. Due to unknown Two-photon absorbance coefficients in porous material, only one-photon case was used in paper. The same script generates lammps executable file and places it to folder /torun.
Files in torun folder can be run via lammps engine and generate modeling result file into foloder outputs. Examples of lammps scripts and result files could be found into corresponding folders.
allscripts.sh can be used for batch execution of all lammps scripts
output_processing.ipynb is an example of data analysis. The most important results of processing are stored in alldata.pkl file (in form of serialized Pandas Dataframe). The way of can be found in output_processing.ipynb
