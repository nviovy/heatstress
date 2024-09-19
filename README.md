# Code for heatstress paper
This repository contains the code used in the paper "Heat stress: an underestimated impact of climate change on vegetation"

1. **ORCHIDEE_heats.tz** contains in a tarball the ORCHIDEE code based on tag 2.2 including the new parameterization of heatstress turnover
2. **CodeC** is a directory than contains the different C codes used to estimate delta Hs_ref, the heatstress indice, from the different remote sensing data (MODIS and COPERNICUS LAI, Ts and different vegetation indices)
3. **R code figures** is a directory that contains both data and R scripts used to generate the figures from the paper

All the codes have been run on linux system on IRENE TGGC supercomputer (https://www-hpc.cea.fr/fr/Joliot-Curie.html)
the ORCHIDEE code is provided with makefiles for compiling the code. So simply launch make on directory ORCHIDEE/src_driver. The others C code can be compiled with standard cc command)
All the codes require the netcdf-4 library and, for ORCHIDEE the OpenMPI library

All the input data to run the codes are too huge (i.e several Tb) to be put on the repository but can be ask on request to nicolas.viovy@lsce.ipsl.fr

