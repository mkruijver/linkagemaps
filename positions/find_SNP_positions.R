
panel_csvs <- list.files("panels", pattern = "\\.csv", full.names = TRUE)

df_panels <- do.call(rbind, lapply(panel_csvs, readr::read_csv))

SNPs <- unique(df_panels$Locus[df_panels$LocusKind == "SNP"])

get_position <- function(snp_name){
  snp_number <- gsub("rs", "", snp_name)

  endpoint <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/", snp_number)

  content <- httr::content(httr::GET(endpoint), encoding = "UTF-8")

  id <- content$present_obs_movements[[1]]$allele_in_cur_release$seq_id

  # dbSNP Variation API positions are 0-based
  # we convert to 1-based GRCh38 coordinates
  pos <- content$present_obs_movements[[1]]$allele_in_cur_release$position + 1L

  data.frame(snp = snp_name, id = id, pos = pos)
}


positions <- list()

i = 1
for (i in seq_along(SNPs)){
  cat(i, "\n")
  positions[[i]] <- get_position(SNPs[i])

  # dbSNP asks not to query the API more than once a second
  Sys.sleep(1)
}

df_positions <- do.call(rbind, positions)
df_positions$chromosome <- as.integer(sub(".*NC_(\\d+).*", "\\1", df_positions$id))

# verify that we pulled the positions corresponding to version 38
grch38 <- readr::read_csv("GRCh38.p14.csv") # from http://togogenome.org/organism/9606
stopifnot(all(df_positions$id %in% grch38$RefSeq))

readr::write_csv(df_positions, "positions/SNP_positions.csv")
