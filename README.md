# README

This repository contains all the model code and plot routines to reproduce the results from Deutloff et al. (2023). 
The repository is structured as follows: 

| Folder        | Description                                                                                                                                                                                                                            |
|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Calibrate     | Contains the calibration routines used to calibrate *r* for PFTP, PFAT and AMAZ. It also contains the scripts to determine the probability distributions used to sample parameters of the model ensembles from.                        |
| Evaluate_Runs | Contains script to calculate aggregated measures of probabilities of triggering and a script to collect data for Table 2.                                                                                                              |
| FaIR          | Contains an adapted version of version 2.0.0 of FaIR which was originally designed by [Leach et al. (2021)](https://gmd.copernicus.org/articles/14/3007/2021/). The adapted version we intruduce here allows for the coupling of CTEM. |
| Plot          | Contains all plotscripts. Which figure a certain plotscript is for is stated in the string at the top of every script.                                                                                                                 |
| Read_data     | Contains routines to read model output data.                                                                                                                                                                                           |
| Runs          | Contains script to run FaIR either coupled or uncoupled to CTEM and produce model ensembles with 5000 members.                                                                                                                         |

To be able to reproduce the results from Deutloff et al. (2024) you also have to download the model output and the 
parameter sets of the model which are provided via Zenodo [here](https://doi.org/10.5281/zenodo.10820467).