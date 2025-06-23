# Wavelength-resolved-magnetic-beam-spin-encoding
Data and source code for "Overcoming velocity broadening effects in magnetic beam encoded microscopy using a wavelength-resolved imaging scheme."

**Please note the analysis code here uses MATLAB version 2022b**

Hi there,

To use the main analysis script ("read_experimental_data.mlx") to analyse the raw experimental data in the study presented, please follow the below instructions;

Download all items inside of this repository and ensure all items are placed into the same directory path. A list of each item is given here;

    "read_experimental_data.mlx" - The main analysis code.
    (FOLDER) "raw datafiles" - Contains 200 .dat files which contain all of the raw experimental data.
    (FOLDER) prereq_echo_reader - Contains a few scripts which are necessary to running 'read_experimental_data.mlx'.
    "phase correction shifts.txt" - Some additional information required to re-phase the experimental data.

On line 15 of "read_experimental_data.mlx" please amend the variable 'pathway' to reflect the directory you've installed all of the above items into.

