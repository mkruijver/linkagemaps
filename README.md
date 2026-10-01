# linkagemaps

R scripts for generating linkage maps for forensic DNA loci.

The maps are based on GRCh38 physical positions and the sex-averaged recombination map from Halldórsson et al. (2019), *Characterizing mutagenic effects of recombination through a sequence-level genetic map*, Science 363:eaau1043. https://doi.org/10.1126/science.aau1043

Currently included:

- **ForenSeq DNA Signature Prep** autosomal STRs and identity SNPs.

## Repository structure

- `generate_maps.R` — generates linkage-map CSV files.
- `panels/` — panel definitions and scripts used to construct them.
- `positions/` — GRCh38 locus positions and scripts used to obtain them.
- `data-raw/aau1043_datas3.gz` — sex-averaged recombination map from Halldórsson et al. (2019).

## Running locally

Run `generate_maps.R` to generate the linkage-map CSV file. Currently only panel is included.

## Licence

Note that the file `data-raw/halldorsson2019/aau1043_datas3.gz` is supplementary data from https://doi.org/10.1126/science.aau1043. The file is included to make generation of the linkage maps reproducible. It is third-party data and is not covered by this repository's MIT licence.