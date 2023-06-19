# panelCAT ReadMe (In progress!)
PanelCAT is an open-source tool to analyse, visualise and compare NGS panel target regions. 
PanelCAT uses R and ShinyR, and is provided per AGPL-3 license, RESEARCH USE ONLY.
PanelCAT is also available via (https://aoszwald.shinyapps.io/panelcat)
## preparation
Clone the repository into an R Project
Download the COSMIC Cancer Mutation Census data from https://cancer.sanger.ac.uk/cosmic/download (cmc_export.tsv) and place it in the db_ori folder
## operation
Run app.R
Go to Tab "NewPanel" -> update All
PanelCat will then acquire the current RefSeq and Clinvar databases, and process the Cosmic CMC database. 
The cmc_export.tsv file can be deleted after it has been processed.
PanelCat will then update all pre-supplied panel analyses.

Panel analyses files can be manually deleted or added from the "panels" folder. The panels are loaded upon running the app.R script. 
