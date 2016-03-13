###
### Load MgDB into namespace
###

.onLoad <- function(libname, pkgname)
{
    ns <- asNamespace(pkgname)
    seq_file <- system.file("extdata", 'rdp_11.4_seq.rds',
                            package=pkgname, lib.loc=libname)

    db_taxa_file <- system.file("extdata", "rdp_11.4.sqlite3",
                                package=pkgname, lib.loc=libname)

    if(!file.exists(seq_file) || !file.exists(db_taxa_file)){
        packageStartupMessage("RDP 11.4 database data not present, use `get_rdpMgDb.R` In the package inst/scripts directory to downlod the database into the package inst/extdata/ directory and reinstall the package")
    }

    metadata = list(URL = "https://rdp.cme.msu.edu",
                    DB_TYPE_NAME = "RDP",
                    DB_VERSION = "rdp_11.4",
                    ACCESSION_DATE = "March 11, 2016")

    ## load database sequence object
    db_seq <- readRDS(seq_file)

    ## initiate new MgDB object
    rdpMgDb <- new("MgDb",
                  seq = db_seq,
                  taxa_file = db_taxa_file,
                  tree_file = "not available",
                  metadata = metadata)

    assign("rdp11.4MgDb", rdpMgDb, envir=ns)
    namespaceExport(ns, "rdp11.4MgDb")

}
