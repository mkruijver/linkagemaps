# linkagemaps

R scripts for generating linkage maps for forensic DNA loci.

The maps are based on GRCh38 physical positions and the sex-averaged recombination map from Halldórsson et al. (2019), *Characterizing mutagenic effects of recombination through a sequence-level genetic map*, Science 363:eaau1043. https://doi.org/10.1126/science.aau1043

Currently included:

- **ForenSeq DNA Signature Prep** autosomal STRs and identity SNPs.

## Data sources and coordinate conventions

The sex-averaged recombination map is taken from the supplementary data of Halldórsson et al. (2019). The file `data-raw/aau1043_datas3.gz` is included in this repository to make generation of the linkage maps reproducible.

Physical positions of SNP loci are obtained from the [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) Variation API. The RefSeq accession returned by dbSNP is checked against the GRCh38 reference sequences in `GRCh38.p14.csv`.

The `position` field returned by the dbSNP API is zero-based. It is therefore converted to a conventional one-based GRCh38 coordinate by adding 1 before the position is stored in `positions/SNP_positions.csv`.

Physical positions of STR loci are obtained from [NIST STRBase](https://strbase.nist.gov/) using the [`STRBaseclient`](https://github.com/mkruijver/STRBaseclient) R package. The script used to retrieve these positions is included as `positions/find_STR_positions.R`.

For STRs, the physical position used is the midpoint of the GRCh38 repeat-region 
coordinates from STRBase. Earlier maps used the STR start instead of the midpoint.

### SNP coordinate correction

The original ForenSeq linkage map released on 14 July 2026 (`v140726`) used the dbSNP API `position` field directly and therefore recorded SNP physical positions one base lower than their corresponding one-based GRCh38 coordinates.

Subsequent maps correct this by adding 1 to the positions returned by the dbSNP API. The correction has only a very small effect on the interpolated genetic positions, but ensures that the physical coordinates use the intended convention. The `v140726` release is retained unchanged so that the historical map remains reproducible.

## Repository structure

- `generate_maps.R` — generates linkage-map CSV files.
- `panels/` — panel definitions and scripts used to construct them.
- `positions/` — GRCh38 locus positions and scripts used to obtain them.
- `data-raw/aau1043_datas3.gz` — sex-averaged recombination map from Halldórsson et al. (2019).

## Running locally

Run `generate_maps.R` to generate the linkage-map CSV file. Currently only panel is included.

## Licence

The code in this repository is released under the MIT licence.

Note that the file `data-raw/aau1043_datas3.gz` is supplementary data from https://doi.org/10.1126/science.aau1043. The file is included to make generation of the linkage maps reproducible. It is third-party data and is not covered by this repository's MIT licence.
