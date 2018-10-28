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
library(R.utils)
library(data.table)
library(digest)

## Database URL
db_root_url <- "https://rdp.cme.msu.edu/download"
seq_bacteria_url <- paste0(db_root_url, "/current_Bacteria_unaligned.fa.gz")
seq_archaea_url <- paste0(db_root_url, "/current_Archaea_unaligned.fa.gz")

## RNAcentral ids - external to RDP
rnacentral_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz"
rnacentral_md5_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/md5/md5.tsv.gz"

## Downloaded files
seq_bacteria_file <- tempfile()
seq_archaea_file <- tempfile()
rnacentralgz_file <- tempfile()
rnacentral_file <- tempfile()
rnacentralgz_md5_file <- tempfile()
rnacentral_md5_file <- tempfile()

## MD5 check sums from initial download
seq_bacteria_md5 <- "c5f43e1285060bde2ca59a26bbd8b824"
seq_archaea_md5 <- "8b93d23c1eb6d96d02e8d53ceef0598e"
rnacentral_md5 <- "41dc9e615933f35e22e2699089872282"

## MgDb database files name
db_file <- "../extdata/rdp11.5.sqlite"
metadata_file <- "../extdata/rdp11.5_metadata.RDS"

### Download database files ####################################################
download_db <- function(url, file_name, md5){
    ## Download file and check to make sure MD5 checksum matches checksum for
    ## previously downloaded version

    download.file(url,file_name)
    new_md5 <- tools::md5sum(file_name)
    if (md5 != new_md5) warning("checksum does not match downloaded file.")
}

## Seq Data
download_db(seq_bacteria_url, seq_bacteria_file, seq_bacteria_md5)
download_db(seq_archaea_url, seq_archaea_file, seq_archaea_md5)

## RNAcentral data
##Directly downloading RNAcentral mapping file
download.file(rnacentral_url, rnacentralgz_file)
gunzip(filename = rnacentralgz_file, destname = rnacentral_file)

##Downloading RNAcentral id to md5 mapping file
download.file(rnacentral_md5_url, rnacentralgz_md5_file)
gunzip(filename = rnacentralgz_md5_file, destname = rnacentral_md5_file)

### Create SQLite DB with Taxa and Seq Data ###################################
### Parse RDP database files - Load seq files,
### seq files contain sequence and taxonomic lineage
parse_rdp <- function(seq){
    rdp_names <- names(seq) %>% str_split("\t")

    rdp_ids <- sapply(rdp_names, function(rdpn) {
        rdpn <- unlist(rdpn)
        ids <- str_split(rdpn[1], pattern = " ", n=2, simplify=TRUE)
        lineage <- str_split(str_replace(rdpn[2], "Lineage=",""),
                             pattern = ";", simplify=TRUE)
        c(ids, lineage)
    })

    get_feature <- function(rr, feature) {
        if(isTRUE(which(rr %in% feature) > 0)) {
            return(rr[which(rr %in% feature) -1])
        }
        ""
    }

    ids <- sapply(rdp_ids, function(rdpn) rdpn[1])
    species <- sapply(rdp_ids, function(rdpn) rdpn[2])
    rootrank <- sapply(rdp_ids, function(rdpn) rdpn[3])
    domain <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "domain"))
    phylum <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "phylum"))
    class <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "class"))
    subclass <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "subclass"))
    ord <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "order"))
    suborder <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "suborder"))
    family <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "family"))
    genus <- sapply(rdp_ids, function(rdpn) get_feature(rdpn, "genus"))

    data.frame(Keys=ids, rootrank=rootrank, domain=domain,
               phylum=phylum, class=class, subclass=subclass,
               ord=ord, suborder=suborder, family=family,
               genus=genus, species=species)
}

seqs <- c(readDNAStringSet(seq_bacteria_file),
         readDNAStringSet(seq_archaea_file))

## formatting lineage into a table
taxa_tbl <- parse_rdp(seqs)

## renaming seqs to match taxonomy table
names(seqs) <-  taxa_tbl$Keys

## Adding RNAcentral and NCBI_tax ids to taxonomy table
add_rnacentral_mapping <- function(rnacentral_md5_file, rnacentral_file, taxa_tbl, seqs){
    md5mapping <- fread(rnacentral_md5_file, sep = "\t", header = FALSE )
    id_map <- fread(rnacentral_file, sep = "\t", header = FALSE )
    new_idmap <- id_map[,c("V1","V4")]
    dedup_idmap <- subset(new_idmap,!duplicated(new_idmap$V1))
    dedup_idmap$digest <- md5mapping$V2[match(dedup_idmap$V1,md5mapping$V1)]
    seqsdigest <- sapply(as.character(seqs), digest, algo = "md5",
                         serialize = F)
    seqsdigest_tbl <- as.data.frame(seqsdigest)
    seqsdigest_tbl$Keys <- names(seqsdigest)
    colnames(seqsdigest_tbl) <- c("md5digest", "Keys")
    seqsdigest_tbl$RNAcentralID <- dedup_idmap$V1[match(seqsdigest_tbl$md5digest,
                                                        dedup_idmap$digest)]
    seqsdigest_tbl$NCBItaxonID <- dedup_idmap$V4[match(seqsdigest_tbl$md5digest,
                                                       dedup_idmap$digest)]
    taxa_tbl$RNAcentralID <- seqsdigest_tbl$RNAcentralID[match(taxa_tbl$Keys,
                                                               seqsdigest_tbl$Keys)]
    taxa_tbl$NCBItaxonID <- seqsdigest_tbl$NCBItaxonID[match(taxa_tbl$Keys,
                                                             seqsdigest_tbl$Keys)]
    ## Return as a data.frame
    data.frame(taxa_tbl)
}

taxa_tbl <- add_rnacentral_mapping(rnacentral_md5_file, rnacentral_file,
                                   taxa_tbl, seqs)

## Creating MgDb formated sqlite database
metagenomeFeatures:::make_mgdb_sqlite(db_name = "rdp11.5",
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
