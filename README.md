# panelCAT ReadMe (In progress!)
PanelCAT is an open-source tool to analyse, visualise and compare NGS panel target regions. 
PanelCAT uses R and ShinyR, and is provided per AGPL-3 license, RESEARCH USE ONLY.
PanelCAT is also available via (https://aoszwald.shinyapps.io/panelcat)
## preparation
- Clone the repository into an R Project.
- Download the COSMIC Cancer Mutation Census data from https://cancer.sanger.ac.uk/cosmic/download (cmc_export.tsv) and place it in the db_ori folder.
- This file needs to be downloaded manually because COSMIC downloads are only available to registered users. Please take note of the COSMIC license conditions for private and institutional users. 
- Run app.R
- Go to Tab "NewPanel" -> UPDATE ALL. PanelCat will acquire current RefSeq and Clinvar databases, and process the COSMIC CMC database. Processing the COSMIC database will take quite some time because the gene names will be updated (and this step is slow).
- PanelCat will then update all pre-supplied panel analyses using the updated databases.
- The cmc_export.tsv file may be deleted after this step is complete.
## operation

## operation - continued
- Analysed panels are stored in the "panels" subfolder.
- Every time panels are analysed or updated, a new file or version of a file is created. These files can be manually removed, replaced, etc.
- PanelCAT will load the most recent version of each panel analysis upon start.
- In case the most recent panel analyses were not all analyzed using the same reference databases, a warning will be issued on startup. This problem can be solved by updating all panels.
- In case panel analyses were not all analyzed using the same version of PanelCat, an error will be issued on startup. This problem can be solved by first updating all panels analyzed with the outdated version using the newer panelcat version, and then adding the newer panels. Alternatively, it is possible to perform new analyses of each panel (this may be easier)
