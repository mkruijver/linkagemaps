# The linkage map from Halldorsson et al. (2019) is used to create
# a linkage map in DBLR format
# The input map is available from https://doi.org/10.1126/science.aau1043
#
# aau1043_datas1.gz based on paternal crossovers
# aau1043_datas2.gz based on maternal crossovers
# aau1043_datas3.gz sex-averaged
#
# we will use aau1043_datas3.gz (the sex-averaged map)
# the map contains the following data:
# Chr (chromosome)
# Begin (start point position of interval in GRCh38 coordinates)
# End (end point position of interval in GRCh38 coordinates)
# cMperMb (recombination rate in interval)
# cM (centiMorgan location of ***end point*** of interval)

# Before we can use this map, we need to find the GRCh38 coordinates
# of all STRs and SNPs contained in the ForenSeq kit
#
# STR positions can be pulled from STRbase
#      (see find_STR_positions.R which creates STR_positions.csv)
# SNP positions can be found in dbSNP
#      (see find_SNP_positions.R which creates SNP_positions.csv)
#

# load and tidy up dataset
decode_map <- readr::read_tsv("data-raw/aau1043_datas3.gz", comment = "#")
names(decode_map)[5] <- "EndcM"

decode_map_split <- split(decode_map, decode_map$Chr)
for (i in seq_along(decode_map_split)){
  # add begin (cM) of intervals to each chromosome
  decode_map_split[[i]]$BegincM <- c(0, decode_map_split[[i]]$EndcM[-nrow(decode_map_split[[i]])])
}

# Two problems
#  Chromosome     Locus       PosBp
#  2              rs876724    114974
#  11             rs2076848   134797652

# rs876724 sits at the start of Chromosome 2 before the start of the decode map
# rs2076848 sits at the end of Chromosome 11 after the end of the decode map
#
# we assume that the recombination rate is approximately 1 cM/Mb
# on the interval outside the coverage and extrapolate
#
# the plots show that this assumption is not unreasonable
# plot(decode_map_split$chr2$Begin/1e6, decode_map_split$chr2$BegincM, type="l",
#      xlim=c(0,2), ylim = c(0,2))
# grid()
#
# plot(decode_map_split$chr11$Begin/1e6, decode_map_split$chr11$BegincM, type="l")
# grid()

# obtain an interpolating function mapping function for each chromosome
extend_bp <- 1e6

decode_map_function_by_chr <- lapply(decode_map_split,
                                     function(chr_df) {
                                       x <- chr_df$Begin
                                       y <- chr_df$BegincM

                                       if (extend_bp > 0){
                                         lhs_x <- chr_df$Begin[1] - extend_bp
                                         lhs_y <- chr_df$BegincM[1] - extend_bp/1e6

                                         rhs_x <-
                                           c(chr_df$End[nrow(chr_df)],
                                             chr_df$End[nrow(chr_df)] + extend_bp)

                                         rhs_y <- c(
                                           chr_df$EndcM[nrow(chr_df)],
                                           chr_df$EndcM[nrow(chr_df)] + extend_bp/1e6)

                                         x <- c(lhs_x, chr_df$Begin, rhs_x)
                                         y <- c(lhs_y, chr_df$BegincM, rhs_y)
                                       }

                                       approxfun(x = x, y = y, method = "linear", rule = 1)
                                     })

# read physical positions of STRs and SNPs
str <- readr::read_csv("positions/STR_positions.csv")
snp <- readr::read_csv("positions/SNP_positions.csv")


panel_name <- "ForenSeq DNA Sig Autosomal STR iSNP"
df_panel <- readr::read_csv(paste0("panels/", panel_name, ".csv"))

str <- str[str$locus %in% df_panel$Locus,]
snp <- snp[snp$snp %in% df_panel$Locus,]

# combine STRs and SNPs in one DataFrame
linkage_map_df <- rbind(data.frame(Chromosome = str$chromosome,
                                   Locus = str$locus,
                                   PosBp = str$start),
                        data.frame(Chromosome = snp$chromosome,
                                   Locus = snp$snp,
                                   PosBp = snp$pos))

# sort by chromosome, then by position
linkage_map_df <- linkage_map_df[order(linkage_map_df$Chromosome, linkage_map_df$PosBp),]

linkage_map_df$PoscM <- sapply(seq_len(nrow(linkage_map_df)), function(i){
  chromosome_number <- paste0("chr", linkage_map_df$Chromosome[i])
  decode_map_function_by_chr[[chromosome_number]](linkage_map_df$PosBp[i])
})


# shift the linkage map by 1 cM for chromosome 2 to avoid a negative cM value
linkage_map_df_shifted <- linkage_map_df
linkage_map_df_shifted$PoscM[linkage_map_df$Chromosome==2] <-
  linkage_map_df$PoscM[linkage_map_df$Chromosome==2] + 1

linkage_map <- data.frame(Chromosome = linkage_map_df_shifted$Chromosome,
                                  Locus = linkage_map_df_shifted$Locus,
                                  `Position (cM)` = linkage_map_df_shifted$PoscM, check.names = FALSE)

# add spaces to Penta D/E locus names
linkage_map$Locus <- gsub(pattern = "PentaE",
                                  replacement = "Penta E", linkage_map$Locus)
linkage_map$Locus <- gsub(pattern = "PentaD",
                                  replacement = "Penta D", linkage_map$Locus)


readr::write_csv(linkage_map,
                 file = "maps/ForenSeq DNA Sig Autosomal STR iSNP_01102026.csv")

