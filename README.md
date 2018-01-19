This is a Bioconductor Annotation package with the Ribosomal Database Project Release 11.5 for use with the metagenomeFeatures package (https://github.com/HCBravoLab/metagenomeFeatures.git).

The package is still in development and not available from bioconductor.

To install the development version of the package:  
1. install metagenomeFeatures, this can be done using devtools `install_github("HCBravoLab/metagenomeFeatures")`  
2. clone this repository `git clone https://github.com/HCBravoLab/ribosomaldatabaseproject11.5MgDb.git`   
3. download the database sequence data and generate the metadata and database file using the `get_rdp11.4MgDb.R` script in `inst/scripts`.

4. install ribosomaldatabaseproject11.4MgDb, this can be done using devtools `install_local("local/path/ribosomaldatabaseproject11.5MgDb")`, replace `local/path` with the path to the downloaded git repo.    

The metagenomeFeatures package has vignettes demonstrating how to work with MgDb-class formated sequence databases.
