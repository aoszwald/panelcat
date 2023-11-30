# PanelCat ReadMe
1. [Disclaimer](#disclaimer)
2. [Introduction](#introduction)
3. [Installation - Windows](#installation---windows)
4. [Installation - Linux](#installation---linux)
5. [Installation on your own Shiny server](#installation-on-your-own-shiny-server)
6. [New panel analysis](#new-panel-analysis)
7. [Advanced](#advanced)
## Disclaimer
- PanelCAT uses specifications of NGS panels (target region files) and public scientific databases (RefSeq, ClinVar, COSMIC) to visualise and compare target regions of NGS Panels. PanelCAT is intended to inform researchers on panel target regions by providing a visualisation of existing data, and to assist in researchers' choice of NGS panels for experimental purposes. The results of PanelCAT are not suggestions for any specific product or manufacturer. The results of this software are for informational use only, research use only, and come with absolutely no guarantee of any form. The results provided by PanelCAT are highly dependent on the input, and on external databases, which are not controlled by PanelCAT.
- The results of this software are NOT intended for medical purposes, e.g., medical advice, diagnosis, prevention, or treatment decisions. As such, this software is NOT intended to assist in a medical setting in determining, e.g., physiological and pathological states, congenital conditions or predisposition, recipient safety or compatibility, treatment response or interactions, or in monitoring therapeutic measures.
- Medical tests may only performed by certified and authorised professionals, using devices in accordance with EU 745/2017 (MDR) or EU 746/2017 (IVDR). This software is intended for research-use only and does NOT constitute a medical device according to EU 745/2017 (MDR), and does NOT constitute an in-vitro diagnostic medial device, or an accessory for an in vitro diagnostic medical device, according to EU 746/2017 (IVDR)
## Introduction
- NGS panels are widely used many areas of biological research. However, correct interpretation of variant calls often requires detailed understanding of the panel target regions. For instance, you may need to know precisely what portion of genes (and more specifically, which exons, or known mutations) are targeted by specific panels. It is useful to directly compare between panels that are available at your lab. PanelCAT (online available at http://panelcat.net) is an open-source tool to analyse, visualise and compare NGS panel target regions. PanelCAT uses R and ShinyR, and is provided per AGPL-3 license. PanelCat is NOT a clinically validated tool and designed for RESEARCH USE ONLY. 
- Regarding BED file format and zero-based indexed coordinates, see https://en.wikipedia.org/wiki/BED_(file_format)#Coordinate_system.
## Installation - Windows
- Install R Statistics (https://cran.rstudio.com/) and R Studio (https://posit.co/download/rstudio-desktop/).
- Perform either of the following (A if you have, or want to use, Git version control software, or B if you do not):
- A) Install Git (https://git-scm.com/download/win), open R Studio and clone the repository into an R Project (File -> New Project -> Version Control -> Git -> Repository URL: https://github.com/aoszwald/panelcat -> Create project as subdirectory of: [Choose a folder you have read and write access to] -> Create Project. 
- B) On the github repository, click the green "code" button, then "Download ZIP". Extract the folder "panelcat_main" and its contents in a location of your choice with read and write access. 
- You may wish to rename the directory from "panelcat_main" to "panelcat" for simplicity.
Open R Studio and create a new R project in the panelcat_main directory that contains App.R: File -> New Project -> Existing Directory -> Browse (and choose the panelcat_main directory containing "App.R" and the rest of the extracted files) -> Create Project.
- Register an account at https://cancer.sanger.ac.uk/cosmic/register.  Please take note of the COSMIC license conditions for private and institutional users. 
- Use your account to log in and download the COSMIC Cancer Mutation Census data from https://cancer.sanger.ac.uk/cosmic/download. Scroll down the list of files to "Cancer Mutation Census" and then press "All data CMC". Be sure to download the Genome GRCh37 version!
- Using your archive software (e.g. 7zip), extract the contents of the downloaded .tar file and then extract the contents of the .gz file contained within. You may need an archive manager tool like 7zip to open the file. In Windows 11, it seems to be possible using the windows folder explorer.
- Place the cmc_export.tsv file into the pre-existent db_ori subdirectory of the PanelCat R project directory. 
- The software has been tested with R version 4.2 and 4.3. If you are using R 4.1 or earlier, you will need to update R (see https://www.r-bloggers.com/2022/01/how-to-install-and-update-r-and-rstudio/ for help). If you are using R 4.2 or 4.3, ensure that your base packages are current by entering the following command in the console window:
```R
update.all(ask = F)
```

- In R Studio, open the file app_full.R. In the right top corner of the window where the opened file is displayed, press the "Run App" button (with the green Arrow). Upon first start, R Studio will install the required packages (This will take a while). Be sure to confirm any prompts about updates with "Yes".
- Once the App starts, it will begin by downloading and processing the RefSeq database. This may take a few minutes. If the network connection is lost during this or the following update steps, the application may terminate unexpectedly.
- Go to Tab "New Analysis" -> UPDATE ALL (at the bottom of the left side panel). PanelCat will acquire current RefSeq and Clinvar databases, and process the COSMIC CMC database. Processing the COSMIC database will take quite some time because the gene names will be updated (and this step is slow).
- PanelCat will then update all pre-supplied panel analyses using the updated databases.
- The cmc_export.tsv file may be deleted after this step is complete.
- The "app_full.R" file contains scripts to update databases and permanently save panel analyses. You should use this file only if you are running it on a local device. However, once you have downloaded current databases, you may also use "app.R" if you do not care about saving new analyses permanently, it will start a little bit faster.
## Installation - Linux
- You can use the same steps as described above for Windows, but you may need to install additional Linux packages first. For example, if using Ubuntu, go into terminal and install the packages libcurl4-openssl-dev, libssl-dev, xml2-config and libxml2-dev, using
```bash 
sudo apt install libcurl4-openssl-dev
sudo apt install libssl-dev
sudo apt install xml2-config
sudo apt install libxml2-dev
```
## Installation on your own Shiny Server
- This is only possible using a Linux environment (you may use a virtual machine, but this is beyond the scope of this guide). Assuming you are using Ubuntu: Download and install shiny server (https://posit.co/download/shiny-server/). 
- The easiest way is to run "app_full.R" on a device, let the app download and process all the databases and panels you wish to provide on your server, and then copy the entire panelcat directory (the one that contains "app.R") into the directory /srv/shiny-server/ on your server.
- By default, shiny server will start one server process listening on port 3838 (this means you will access the panelcat app running on your server by entering e.g. "123.456.789.12:3838/panelcat" (123.456.789.12 being a placeholder for your IP in the network).
- Accessing this adress will activate the script in the file "app.R" (not "app_full.R").
- In contrast to the "app_full.R" script, the "app.R" file does not include functions to update the databases and save new analyses permanently. This is useful on a server, where you want to have full control over updates and permanently saved panels (and not the users who access your server via browser!).
- You *could* delete "app.R", then rename "app_full.R" to "app.R", and then give the shiny app ownership over its own folder: 
```bash
sudo chown shiny:shiny /srv/shiny-server/panelcat
```
- This *may* enable clients to access your server with the full functionality of permanently saving panels, and performing database updates (this has not been tested). *However*, do so at your own risk, since it may compromise the security of your system.
## New panel analysis
- PanelCat uses hg19 (GRCh37) based databases. Please ensure that your target region file is also hg19 (GRCh37) based.
- To analyse a new panel, go to the "NewPanel" tab. Target regions need to be provided as a tab-separated table (typically .bed or .txt) that includes the following columns: Chromosome, Start, Stop (e.g. chr1	27100287	27100394). Additional columns may be present.
- Note that the chromosome ID needs the "chr" prefix, e.g. "chr1" (not "1"), or "chrX" (not "X").
- If the software does not automatically detect the first row of the coordinates, a start row (and optionally, a stop row) may be specified by entering the row numbers seperated by a comma (eg. "5" or "5,250")
- If the first three columns of the table do not exactly match the columns Chromosome, Start, Stop (in this order), specify the column number and order (e.g. "2,3,1")
- Optionally, provide a mask file (again .bed or .txt) in the same manner as the target regions file.
- Press "start". Once completed, the panel analyses will be either permanently stored within the "panels" subfolder (when using "App_full.R"), or only temporarily available within the current R session / Shiny server instance (App.R). When using App.R, you can download the analysis result, and upload it to the software at a later time.
## Advanced
- Analysed panels are stored in the "panels" subfolder.
- Every time panels are analysed or updated, a new file or version of a file is created. These files can be manually removed, replaced, etc.
- PanelCAT will load the most recent version of each panel analysis upon start.
- In case the most recent panel analyses were not all analyzed using the same reference databases, a warning will be issued on startup. This problem can be solved by updating all panels.
- In case panel analyses were not all analyzed using the same version of PanelCat, an error will be issued on startup. This problem can be solved by first updating all panels analyzed with the outdated version using the newer panelcat version, and then adding the newer panels. Alternatively, it is possible to perform new analyses of each panel (this may be easier)
