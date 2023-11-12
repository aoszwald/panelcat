# panelCAT ReadMe (In progress!)
PanelCAT is an open-source tool to analyse, visualise and compare NGS panel target regions. 
PanelCAT uses R and ShinyR, and is provided per AGPL-3 license, RESEARCH USE ONLY.
PanelCAT is also available via (http://panelcat.net)
## Installation - Windows
- Install R Statistics (https://cran.rstudio.com/) and R Studio (https://posit.co/download/rstudio-desktop/)
- Perform either of the following (A if you have, or want to use, Git version control software, or B if you do not):
- A) Install Git (https://git-scm.com/download/win), open R Studio and clone the repository into an R Project (File -> New Project -> Version Control -> Git -> Repository URL: https://github.com/aoszwald/panelcat -> Create project as subdirectory of: [Choose a folder you have read and write access to] -> Create Project. 
- B) On the github repository, click the green "code" button, then "Download ZIP". Extract the folder "panelcat_main" and its contents in a location of your choice with read and write access. Open R Studio and create a new R project in the panelcat_main directory that contains App.R: File -> New Project -> Existing Directory -> Browse (and choose the panelcat_main directory containing "App.R" and the rest of the extracted files) -> Create Project.
- Register an account at https://cancer.sanger.ac.uk/cosmic/register.  Please take note of the COSMIC license conditions for private and institutional users. 
- Use your account to log in and download the COSMIC Cancer Mutation Census data from https://cancer.sanger.ac.uk/cosmic/download. Scroll down the list of files to "Cancer Mutation Census" and then press "All data CMC". Be sure to download the Genome GRCh37 version!
- Using your archive software (e.g. 7zip), extract the contents of the downloaded .tar file and then extract the contents of the .gz file contained within. You may need an archive manager tool like 7zip to open the file. In Windows 11, it seems to be possible using the windows folder explorer.
- Place the cmc_export.tsv file into the pre-existent db_ori subdirectory of the PanelCat R project directory. 
- The software has been tested with R version 4.2 and 4.3. If you are using R 4.1 or earlier, you will need to update R (see https://www.r-bloggers.com/2022/01/how-to-install-and-update-r-and-rstudio/ for help). If you are using R 4.2 or 4.3, ensure that your base packages are current by entering "update.all(ask = F)" and pressing Enter into the console window.
- In R Studio, open the file app_full.R. In the right top corner of the window where the opened file is displayed, press the "Run App" button (with the green Arrow). Upon first start, R Studio will install the required packages (This will take a while). Be sure to confirm any prompts about updates with "Yes".
- Once the App starts, it will begin by downloading and processing the RefSeq database. This may take a few minutes. If the network connection is lost during this or the following update steps, the application may terminate unexpectedly.
- Go to Tab "New Analysis" -> UPDATE ALL (at the bottom of the left side panel). PanelCat will acquire current RefSeq and Clinvar databases, and process the COSMIC CMC database. Processing the COSMIC database will take quite some time because the gene names will be updated (and this step is slow).
- PanelCat will then update all pre-supplied panel analyses using the updated databases.
- The cmc_export.tsv file may be deleted after this step is complete.
## operation - basic
- Run "App_full.R" (for full functionality) or "App.R" (in case of hosting the app in a network, see details below)
- View analyses results in different tabs. Most of the time, the presentation of results can be controlled by adjusting one or more input parameters.
## operation - New panel analysis
- PanelCat uses hg19 (GRCh37) based databases. Please ensure that your target region file is also hg19 (GRCh37) based.
- To analyse a new panel, go to the "NewPanel" tab. Target regions need to be provided as a tab-separated table (typically .bed or .txt) that includes the following columns: Chromosome, Start, Stop (e.g. chr1	27100287	27100394). Additional columns may be present.
- Note that the chromosome ID needs the "chr" prefix, e.g. "chr1" (not "1"), or "chrX" (not "X").
- If the software does not automatically detect the first row of the coordinates, a start row (and optionally, a stop row) may be specified by entering the row numbers seperated by a comma (eg. "5" or "5,250")
- If the first three columns of the table do not exactly match the columns Chromosome, Start, Stop (in this order), specify the column number and order (e.g. "2,3,1")
- Optionally, provide a mask file (again .bed or .txt) in the same manner as the target regions file.
- Press "start". Once completed, the panel analyses will be either permanently stored within the "panels" subfolder (when using "App_full.R"), or only temporarily available within the current R session / Shiny server instance (App.R). When using App.R, you can download the analysis result, and upload it to the software at a later time.
## operation - advanced
- Analysed panels are stored in the "panels" subfolder.
- Every time panels are analysed or updated, a new file or version of a file is created. These files can be manually removed, replaced, etc.
- PanelCAT will load the most recent version of each panel analysis upon start.
- In case the most recent panel analyses were not all analyzed using the same reference databases, a warning will be issued on startup. This problem can be solved by updating all panels.
- In case panel analyses were not all analyzed using the same version of PanelCat, an error will be issued on startup. This problem can be solved by first updating all panels analyzed with the outdated version using the newer panelcat version, and then adding the newer panels. Alternatively, it is possible to perform new analyses of each panel (this may be easier)

- "App.R" does not include functions to permanently save panel analyses or update panels. This version of the software should be run when using a hosting platform such as shinyserver (https://posit.co/products/open-source/shinyserver/). This is also the version available at (https://aoszwald.shinyapps.io/panelcat)
