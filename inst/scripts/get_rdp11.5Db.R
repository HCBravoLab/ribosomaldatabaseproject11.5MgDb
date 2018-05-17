###############################################################################
## Code to generate ribosomal database project annotation package extdata the
## species level id includes a mix of species and sub-species names, e.g.
## clones, strains, ect.
## R package with code to parse different database file formats
library(tidyverse)
library(stringr)
library(DECIPHER)
library(Biostrings)
library(metagenomeFeatures) ## Needed for the make_mgdb_sqlite function

## Database URL
db_root_url <- "https://rdp.cme.msu.edu/download"
seq_bacteria_url <- paste0(db_root_url, "/current_Bacteria_unaligned.fa.gz")
seq_archaea_url <- paste0(db_root_url, "/current_Archaea_unaligned.fa.gz")

## RNAcentral ids - external to RDP
rnacentral_url <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral",
                         "/releases/8.0/id_mapping/database_mappings/",
                         "rdp.tsv")

## Downloaded files
seq_bacteria_file <- tempfile()
seq_archaea_file <- tempfile()
rnacentral_file <- tempfile()

## MD5 check sums from initial download
seq_bacteria_md5 <- "c5f43e1285060bde2ca59a26bbd8b824"
seq_archaea_md5 <- "8b93d23c1eb6d96d02e8d53ceef0598e"
rnacentral_md5 <- "41dc9e615933f35e22e2699089872282"

## MgDb database files name
db_file <- "../extdata/rdp11.5.sqlite"
metadata_file <- "../extdata/rdp11.5_metadata.RDS"

### Download database files ####################################################
download_db <- function(url, file_name, md5){
    ## Downloade file and check to make sure MD5 checksum matches checksum for
    ## previously downloaded version

    download.file(url,file_name)
    new_md5 <- tools::md5sum(file_name)
    if (md5 != new_md5) warning("checksum does not match downloaded file.")
}

## Seq Data
download_db(seq_bacteria_url, seq_bacteria_file, seq_bacteria_md5)
download_db(seq_archaea_url, seq_archaea_file, seq_archaea_md5)

## RNAcentral data
download_db(rnacentral_url, rnacentral_file, rnacentral_md5)

### Create SQLite DB with Taxa and Seq Data ###################################
### Parse RDP database files
get_rdp_lineage <- function(rdp_names, rdp_ids){
    rdp_lineage <- rdp_names %>% map(2) %>%
        str_replace("Lineage=","") %>% str_split(";")

    name_lineage <- function(lineage){
        lin_rank <- lineage[c(FALSE,TRUE)]
        lin_taxa <- lineage[c(TRUE,FALSE)]
        data_frame(rank = lin_rank, taxa = lin_taxa)
    }

    rdp_lineage %>%
        set_names(flatten_chr(rdp_ids$ids)) %>%
        map_df(name_lineage, .id = "Key") %>%
        spread(rank, taxa)
        mutate(species = rdp_ids$species %>% flatten_chr())
}

get_rdp_lineage <- function(rdp_names, rdp_ids){
    rdp_lineage <- rdp_names %>% map(2) %>%
        str_replace("Lineage=","") %>% str_split(";")

    name_lineage <- function(lineage){
        lin_rank <- lineage[c(FALSE,TRUE)]
        lin_taxa <- lineage[c(TRUE,FALSE)]
        data_frame(rank = lin_rank, taxa = lin_taxa) %>%
            filter(rank != "") %>% spread(rank, taxa)
    }

    rdp_lineage %>% map_df(name_lineage, .id = "Key") %>%
    mutate(Key = rdp_ids$ids %>% flatten_chr(),
           species = rdp_ids$species %>% flatten_chr())
}

get_rdp_ids <- function(rdp_names){
    rdp_names %>% map(1) %>%
        str_split(pattern = " ", n = 2) %>%
        transpose() %>% set_names(c("ids","species"))
}

### Load seq files, seq files contain sequence and taxonomic lineage
parse_rdp <- function(seq){
    rdp_names <- names(seq) %>% str_split("\t")

    rdp_ids <- rdp_names %>%
        map(1) %>%
        str_split(pattern = " ", n = 2) %>%
        transpose() %>%
        set_names(c("ids","species")) %>%
        map(flatten_chr)

    rdp_lineage <- rdp_names %>%
        map(2) %>%
        str_replace("Lineage=","") %>%
        str_split(";")

    name_lineage <- function(lineage){
        data.frame(rank = lineage[c(FALSE,TRUE)],
                   taxa = lineage[c(TRUE,FALSE)],
                   stringsAsFactors = FALSE)
    }

    lineage_df <- rdp_lineage %>%
        set_names(rdp_ids$ids) %>%
        map_df(name_lineage, .id = "Keys")

    ## Returns wide data frame with desired taxa ordering
    lineage_df %>%
        mutate(taxa = str_trim(taxa, side = "both"),
               taxa = str_replace_all(taxa, '\\"', "")) %>%
        spread(rank, taxa) %>%
        mutate(species = rdp_ids$species) %>%
        select(Keys, rootrank, domain, phylum, class, subclass,
               order, suborder, family, genus, species) %>%
        ## Rename due to SQLite funciton name conflict
        dplyr::rename(ord = order)
}

seqs <- c(readDNAStringSet(seq_bacteria_file),
         readDNAStringSet(seq_archaea_file))

## formatting lineage into a table
taxa_tbl <- parse_rdp(seqs)

## creating seq RDS
names(seqs) <-  taxa_tbl$Keys

## Load RNAcentral data
rnacentral_ids <- read.delim(rnacentral_file,
                             stringsAsFactors = FALSE,
                             header = FALSE)

colnames(rnacentral_ids) <- c("rnacentral_ids", "database",
                              "Keys", "ncbi_tax_id",
                              "RNA_type", "gene_name")

## Dropping database, RNA_type, and gene_name columns
rnacentral_ids$database <- NULL
rnacentral_ids$RNA_type <- NULL
rnacentral_ids$gene_name <- NULL

## Converting ids to character strings - ensure compatible typing
rnacentral_ids$Keys <- as.character(rnacentral_ids$Keys)
taxa_tbl$Keys <- as.character(taxa_tbl$Keys)

## Adding RNAcentral and NCBI_tax ids to taxonomy table
taxa_tbl <- dplyr::left_join(taxa_tbl, rnacentral_ids)

## Creating MgDb formated sqlite database
metagenomeFeatures::make_mgdb_sqlite(db_name = "rdp11.5",
                                     db_file = db_file,
                                     taxa_tbl = taxa_tbl,
                                     seqs = seqs)

### Database Metadata #########################################################
rdp_metadata <- list(ACCESSION_DATE = date(),
                 URL = "https://rdp.cme.msu.edu",
                 DB_TYPE_NAME = "RDP",
                 DB_VERSION = "11.5",
                 DB_TYPE_VALUE = "MgDb",
                 DB_SCHEMA_VERSION = "2.0")

saveRDS(rdp_metadata, file = metadata_file)
