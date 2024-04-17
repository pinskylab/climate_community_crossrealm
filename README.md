# climate_community_crossrealm
Analyses of the relationship between temporal turnover in species composition and temperature change across realms.

Supporting data and code for:  
*Warming and cooling catalyze widespread temporal turnover in biodiversity*  
by Malin L. Pinsky, Helmut Hillebrand, Shane A. Blowes, Jonathan M. Chase, Laura H. Antão, Myriam R. Hirt, Ulrich Brose, Michael T. Burrows, Benoit Gauzens, and Benjamin Rosenbaum


## Basic repo structure
*see readmes in each directory for further information*
- code/: has scripts for data analysis and producing figures and tables
- data/: has data unique to this project
- dataDL/: has data downloaded from other sources. Not tracked by Git
  - bowler_atcs: Human impact data from [Bowler et al. 2020](https://doi.org/10.1002/pan3.10071). Specifically the unzipped contents of pan310071-sup-0003-Supinfo2.7z from the SOM.
  - ersst/: ERSST v5 from ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
  - cruts/:  CRU TS v.4.03 from http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.03/data/cru_ts4.03.1901.2018.tmp.dat.nc
  - worldclim/: WorldClim climatology from https://worldclim.org (wc2.0_bio_30s_01.tif)
- figures/: plots produced by scripts
- output/: files produced by scripts
- temp/: temporary files created by scripts from data in data/ or dataDL/. Not tracked by Git
