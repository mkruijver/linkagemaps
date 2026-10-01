panel_csvs <- list.files("panels", pattern = "\\.csv", full.names = TRUE)

df_panels <- do.call(rbind, lapply(panel_csvs, readr::read_csv))

STRs <- unique(df_panels$Locus[df_panels$LocusKind == "STR"])

require(STRBaseclient)

positions <- list()

for (i in seq_along(STRs)){
  cat(i, "\n")

  locus <- STRs[i]
  lookup <- STRBaseclient::get_locus_chromosome_location_and_allele_reference_table(locus)

  i_row <- grep("38", lookup$assembly)

  positions[[i]] <- lookup[i_row,]

  Sys.sleep(0.5)
}

df_positions <- do.call(rbind, positions)

readr::write_csv(df_positions, "positions/STR_positions.csv")
