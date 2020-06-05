# climate_community_crossrealm
Comparisons of community responses to temperature change across realms

Basic repo structure
- code/: has scripts for data analysis
- data/: has data unique to this project
- dataDL/: has data downloaded from other sources. Not tracked by Git
  - ersst/: ERSST v5 from ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
  - cruts/:  CRU TS v.4.03 from http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.03/data/cru_ts4.03.1901.2018.tmp.dat.nc
  - eltontraits/: Elton Traits from http://www.esapubs.org/archive/ecol/E095/178/#data
  - ledatraits/: LEDA traits from https://uol.de/en/landeco/research/leda/data-files/
  - ocean_productivity/: Ocean productivity from http://orca.science.oregonstate.edu/2160.by.4320.yearly.hdf.land.ocean.merge.php
  - try/: TRY traits from a data request to https://www.try-db.org
  - worldclim/: WorldClim climatology from https://worldclim.org (wc2.0_bio_30s_01.tif)
- output/: files produced by scripts
- figures/: plots produced by scripts
- temp/: temporary files created by scripts from data in data/ or dataDL/. Not tracked by Git
