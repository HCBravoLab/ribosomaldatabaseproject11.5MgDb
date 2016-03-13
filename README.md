This is an R package with the Ribosomal Database Project Release 11.4 for annotating 16S metagenomic sequence data using the metagenomeFeatures package (https://github.com/HCBravoLab/metagenomeFeatures.git).

The package is still in development and not available from bioconductor.

To install the development version of the package:  
1. install metagenomeFeatures, this can be done using devtools `install_github("HCBravoLab/metagenomeFeatures")`  
2. clone this repository `git clone https://github.com/HCBravoLab/ribosomaldatabaseproject11.4MgDb.git`   
3. to download the database data, from the `inst/scripts` directory in the repository run the `get_ribosomaldatabaseproject11.4MgDb.R` script. Note the script requires `dplyr`, `RSQLite`, `phyloseq`, and `Biostrings`.  

4. install ribosomaldatabaseproject11.4MgDb, this can be done using devtools `install_local("local/path/ribosomaldatabaseproject11.4MgDb")`, replace `local/path` with the path to the downloaded git repo.    

The metagenomeFeatures package has vignettes demonstrating how to use the ribosomaldatabaseproject11.4MgDb package to annotate 16S rRNA metagenomic sequence data.
