From Michael Burrows, 29 March 2020

Here’s that taxonomic breakdown of the Biotime studies.

1. biotime_taxa.csv   = all the biotime taxa matched with phylum, class and family where possible.

2. TaxGroupphylumclassrealmtab.csv = Phylum and class with frequency of records in freshwater, marine and terrestrial studies, and a manually added field representing major groups for splitting analyses.

3. TaxGroupByStudy.csv = Study detail with dominant taxon for each study and 2degree lat/long gridcell combination. Proportions of records in each major group are in the trailing columns (7 onwards). Dominant group is that with the highest proportion of records in the study.

This only works for the species-level taxa so far, but those do make up most of the records.


From Michael Burrows, 21 May 2002

1. stibysprealm.csv: Here is my STI data for all the BIOTIME species, including NAs for species that failed to fit Maxent models to OBIS/GBIF data, usually because of insufficient data. The percentile temperatures (X0, X10, … X100) are derived using Maxent predicted global distributions to extract data from average surface temperatures for land and ocean.