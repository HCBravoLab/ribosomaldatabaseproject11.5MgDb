###
### Load MgDB into namespace
###
.onAttach <- function(libname, pkgname){

    db_file <- system.file("extdata", "rdp11.5.sqlite",
                           package = pkgname, lib.loc = libname)

    metadata_file <- system.file("extdata", "rdp11.5_metadata.RDS",
                                 package = pkgname, lib.loc = libname)

    ## Note no tree for rdp11.5

    if (!file.exists(db_file) | !file.exists(metadata_file)) {
        packageStartupMessage("RDP 11.5 database data not present,
                              use `get_rdp11.5Db` In the package inst/scripts
                              directory to download the database into the
                              package inst/extdata/ directory and reinstall
                              the package")
    }
}

.onLoad <- function(libname, pkgname){
    ns <- asNamespace(pkgname)

    db_file <- system.file("extdata", "rdp11.5.sqlite",
                           package = pkgname, lib.loc = libname)

    metadata_file <- system.file("extdata", "rdp11.5_metadata.RDS",
                                 package = pkgname, lib.loc = libname)

    ## Note no tree for rdp11.4

    metadata <- readRDS(metadata_file)

    ## initiate new MgDB object
    rdpMgDb <- metagenomeFeatures::newMgDb(db_file = db_file,
                                          tree = NULL,
                                          metadata = metadata)

    assign("rdp11.5MgDb", rdpMgDb, envir = ns)
    namespaceExport(ns, "rdp11.5MgDb")

}

