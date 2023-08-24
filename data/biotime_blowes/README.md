# climate_community_crossrealm/data
Contains files processed from BioTime.

- all_pairs_beta.Rdata: pairwise dissimilarity between each observation within each BioTime timeseries after gridding to 96 km<sup>2</sup> and removal of observations with <85% coverage. Has Jaccard replacement component dissimilarity (Jtu) and Morisita-Horn similarity (Horn).
- bt_grid_spp_list_abund.Rdata: Has two objects (bt_grid_spp_list_abund and bt_grid_spp_list_abund_year). The first has total abundance (all samples, all years). The second has total annual abundance (all samples).
- bt_grid_spp_list.Rdata: species list for each time series
- bt_malin.Rdata: has metadata on each timeseries, including
  - REALM: one of terrestrial, freshwater or marine
  - Biome: the biomes used in the Blowes et al. 2019 Science analysis
  - taxa_mod: the taxonomic groups used in the 2019 analysis (and are somewhat simplified from the BioTIME database)
  - STUDY_ID: an identifier from BioTIME
  - rarefyID: concatenation of STUDY_ID and the cell reference number from the gridding process
  - rarefyID_x: the longitude of the centroid from a convex hull around the data within cells
  - rarefyID_y: the latitude of the centroid from a convex hull around the data within cells
  - S is species richness
- time_series_data_type.Rdata: has rarefyID_type with two columns
  - rarefyID: see bt_malin.Rdata
  - BROAD_TYPE: one of "count", "presence", or "biomass" to describe the type of data in the timeseries
