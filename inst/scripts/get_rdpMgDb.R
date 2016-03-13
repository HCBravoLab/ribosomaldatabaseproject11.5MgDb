library(Biostrings)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)
# library(tidyr)

name_lineage <- function(lineage){
    lin_rank <- lineage[c(FALSE,TRUE)] %>% paste0("r_", .)
    lineage[c(TRUE,FALSE)] %>%
        as.list() %>% set_names(lin_rank) %>%
        as_data_frame()
}

get_rdp_lineage <- function(rdp_names, rdp_ids){
    rdp_lineage <- rdp_names %>% map(2) %>%
        str_replace("Lineage=","") %>% str_split(";")
    rdp_lineage %>% map_df(name_lineage) %>% mutate(Key = rdp_ids)
}

get_rdp_ids <- function(rdp_names){
    rdp_names %>% map(1) %>%
        str_split(pattern = " ", n = 2) %>%
        map(1) %>% flatten_chr()
}

create_rdp_db <- function(fasta_file, db_name = "rdp_11.4"){
    seq <- readDNAStringSet(fasta_file)

    ## extracting RDP ids
    rdp_names <- names(seq) %>% str_split("\t")
    rdp_ids <- get_rdp_ids(rdp_names)

    ## formatting lineage into a table
    taxa <- get_rdp_lineage(rdp_names, rdp_ids)

    ## creating seq RDS
    names(seq) <- rdp_ids
    db_seq_file <- paste0("../extdata/",db_name, "_seq.rds")
    saveRDS(seq,db_seq_file)

    ## creating taxa sqlite
    db_taxa_file <- paste0("../extdata/",db_name, ".sqlite3")
    db_con <- dplyr::src_sqlite(db_taxa_file, create = T)
    dplyr::copy_to(db_con,taxa,
                   temporary=FALSE,
                   indexes=list(colnames(taxa)))

}

## current_Prokayote_unaligned_split*fa.gz includes Both Archaea and Bacteria 16S rRNA sequences from RDP release 11.4
##
## https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
system("cat ~/Downloads/current_Archaea_unaligned.fa ~/Downloads/current_Bacteria_unaligned.fa > ~/Desktop/current_Prokaryote_unaligned.fa")
fasta_file <- "current_Prokaryote_unaligned.fa.gz"
create_rdp_db(fasta_file)
