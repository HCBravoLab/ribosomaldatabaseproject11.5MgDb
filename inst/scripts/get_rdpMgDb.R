## Code to generate ribosomal database project annotation package extdata the
## species level id includes a mix of species and sub-species names, e.g.
## clones, strains, ect.

library(Biostrings)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)
library(tidyr)

name_lineage <- function(lineage){
    lin_rank <- lineage[c(FALSE,TRUE)]
    lin_taxa <- lineage[c(TRUE,FALSE)]
    data_frame(rank = lin_rank, taxa = lin_taxa) %>%
        filter(rank != "") %>% spread(rank, taxa)
}

get_rdp_lineage <- function(rdp_names, rdp_ids){
    rdp_lineage <- rdp_names %>% map(2) %>%
        str_replace("Lineage=","") %>% str_split(";")

    rdp_lineage %>% map_df(name_lineage) %>%
        mutate(Key = rdp_ids$ids %>% flatten_chr(),
               species = rdp_ids$species %>% flatten_chr())
    ## code to parse species names
    #%>% separate(species,into = c("species","clone"), sep = ";")
}

get_rdp_ids <- function(rdp_names){
    rdp_names %>% map(1) %>%
        str_split(pattern = " ", n = 2) %>%
        transpose() %>% set_names(c("ids","species"))
}

create_rdp_db <- function(fasta_file, db_name = "rdp_11.4"){
    seq <- readDNAStringSet(fasta_file)

    ## extracting RDP ids
    rdp_names <- names(seq) %>% str_split("\t")
    rdp_ids <- get_rdp_ids(rdp_names)

    ## formatting lineage into a table
    taxa <- get_rdp_lineage(rdp_names, rdp_ids)

    ## creating seq RDS
    names(seq) <- rdp_ids$ids %>% flatten_chr()
    db_seq_file <- paste0("../extdata/",db_name, "_seq.rds")
    saveRDS(seq,db_seq_file)

    ## creating taxa sqlite
    db_taxa_file <- paste0("../extdata/",db_name, ".sqlite3")
    db_con <- dplyr::src_sqlite(db_taxa_file, create = T)
    dplyr::copy_to(db_con,taxa,
                   temporary=FALSE,
                   indexes=list(colnames(taxa)))

}

## current_Prokayote_unalignedfa.gz includes Both Archaea and Bacteria 16S rRNA
## sequences from RDP release 11.4 use download_rdp.sh to dowload and
## concatenate the Bacteria and Archaea RDP seq files
fasta_files <- c("current_Prokaryote_unaligned.fa.gz")
create_rdp_db(fasta_files)
