# Dependencies ------------------------------------------------------------
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)}
if(!require(shinyjs)){
  install.packages("shinyjs")
  library(shinyjs)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)}
if(!require(plotly)){
  install.packages("plotly")
  library(plotly)}
if(!require(DT)){
  install.packages("DT")
  library(DT)}
if(!require(shinyalert)){
  install.packages("shinyalert")
  library(shinyalert)}
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)}
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)}
if(!require(randomcoloR)){
  install.packages("randomcoloR")
  library(randomcoloR)}
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.17")} 
if(!require(GenomicFeatures)){
  BiocManager::install("GenomicFeatures")
  library(GenomicFeatures)}
if(!require(fst)){
  install.packages("fst")
  library(fst)}
if(!require(ggrepel)){
  install.packages("ggrepel")
  library(ggrepel)}
options(repos = BiocManager::repositories())
options(timeout = 360)

# Functions ---------------------------------------------------------------
{
  # set true to always skip database update
  debugMode <- F
  
  # PanCatversion
  PanCatv <- 24
  
  count_nna_func <- function(x) sum(!is.na(x))
  
  # reverse plotly legend labels
  reverse_legend_labels <- function(plotly_plot) {
    n_labels <- length(plotly_plot$x$data)
    plotly_plot$x$data[1:n_labels] <- plotly_plot$x$data[n_labels:1]
    plotly_plot
  }
  
  # for loading dbx names
  nameList <- function(lst, n){
    sapply(lst, `[`, n)
  }
  
  # for processing clinVar
  remove_gene_number <- function(x) {
    str_split_fixed(x, ":",2)[,1]
  }
  
  # RefSeq loader
  loadRefSeq <- function(updateDb) {
    if ((updateDb == FALSE | debugMode == TRUE) & length(txdb_path) != 0) {
      if (exists("txdb") == F) {
        incProgress(1, detail = "Loading RefSeq database")
        txdb <<- loadDb(txdb_path)
      } else {
        incProgress(1, detail = "Loading RefSeq database")
      }
    } else {
      incProgress(1, detail = "update RefSeq db")
      txdb <<- makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz",
                               dataSource = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz")
      incProgress(0, detail = "saving txdb")
      saveDb(txdb, paste0("txdb_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
      txdb_path <<- list.files(pattern = "txdb_", full.names = T)[length(list.files(pattern = "txdb_"))]
    }
    if (length(assRep_path != 0)) {
      assRep <<- readRDS(assRep_path)
    } else {
      assRep <<- read.table("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt", fill = T)
      saveRDS(assRep, paste0("assRep_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
    }
    seqlevels(txdb) <<- seqlevels(txdb)[1:24]
    seqlevels(txdb) <<- assRep$V11[1:24] 
    if (exists("ex_by_ge") == F) {
      incProgress(0, detail = "Get Exons")
      ex_by_ge <<- exonsBy(txdb, by="gene")
    }
  }
  
  createExonTable <- function() {
    if (!exists("ex_by_ge_df1")) {
      withProgress(message = "Preparation", value = 0, max = 4, {
        loadRefSeq(updateDb = FALSE)
        ex_by_ge_df1 <- (as.data.table(unlist(ex_by_ge)))
        splitnames <- str_split_fixed(ex_by_ge_df1$exon_name, "-",3)
        colnames(splitnames) <- c("name","transcript","exon")
        ex_by_ge_df1 <- cbind(ex_by_ge_df1, splitnames)
        ex_by_ge_df1 <- subset(ex_by_ge_df1, select = -c(exon_name,name))
        ex_by_ge_df1$exon <- as.numeric(ex_by_ge_df1$exon)
        ex_by_ge_df1 <<- ex_by_ge_df1
      })
    }
  }
  
  # CLINVAR update check and loader
  updateCheckClv <- function() {
    url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz.md5"
    tmp <- tempfile()
    download.file(url, tmp)
    clvMd5Serv <- readLines(tmp)
    
    if (length(clv_md5_path) == 0 ){
      return(TRUE)
    } else {
      clvMd5Cli <- readRDS(clv_md5_path)
      if (identical(clvMd5Serv, clvMd5Cli) == T) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
  }
  
  loadClinVar <- function(updateDb) {
    if ((updateDb == FALSE | debugMode == TRUE) & length(clv_path) != 0) {
      if (exists("clinvar") == F) {
        incProgress(1, detail = "Loading ClinVar database")
        clinvar <<- read.fst(clv_path, as.data.table = T)
      } else {
        incProgress(1, detail = "Loading ClinVar database")
      }
    } else {
      incProgress(1, detail = "Updating ClinVar db")
      url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz.md5"
      tmp <- tempfile()
      download.file(url, tmp)
      clvMd5Serv <- readLines(tmp)
      saveRDS(clvMd5Serv, paste0("clvmd5_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
      
      url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz"
      tmp <- tempfile()
      download.file(url, tmp)
      
      clinvar <- as.data.table(read.table(gzfile(tmp)))
      colnames(clinvar) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
      clinvar$QUAL <- NULL
      clinvar$FILTER <- NULL
      clinvar$`#CHROM` <- paste0("chr",clinvar$`#CHROM`)
      clinvar$clnsig <- str_match(clinvar$INFO, "CLNSIG=\\s*(.*?)\\s*;")[,2]
      clinvar$clnrevstat <- str_match(clinvar$INFO, "CLNREVSTAT=\\s*(.*?)\\s*;")[,2]
      clinvar <- clinvar[!is.na(clinvar$clnsig),]
      clinvar$gene <- str_match(clinvar$INFO, "GENEINFO=\\s*(.*?)\\s*;")[,2]
      clinvar$INFO <- NULL
      clinvar <- clinvar %>% mutate(gene = str_split(clinvar$gene, "\\|")) %>%
        unnest(gene)
      clinvar$gene <- remove_gene_number(clinvar$gene)
      clinvar <- clinvar[,c("gene","#CHROM","POS","ID","REF","ALT","clnsig","clnrevstat")]
      incProgress(0, detail = "saving clinvar db")
      write.fst(clinvar, paste0("clinvar_", format(Sys.time(), "%Y%m%d_%H%M%S"),".fst"), compress = 100)
      clinvar <<- clinvar
      clv_path <<- list.files(pattern = "clinvar_", full.names = T)[length(list.files(pattern = "clinvar_"))]
    }
  }
  
  # COSMIC loader
  loadCosmic <- function(updateDb) {
    if ((updateDb == FALSE | debugMode == TRUE) & length(cmc_path) != 0) {
      if (exists("cmc_ori") == F) {
        incProgress(1, detail = "reading CMC db")
        cmc_ori <<- read.fst(cmc_path, as.data.table = T)
      } else {
        incProgress(1, detail = "reading CMC db")
      }
    } else {
      library(Seurat)
      incProgress(1, detail = "Get COSMIC db")
      cmc_ori <- fread("db_ori/cmc_export.tsv")
      cmc_genes <- data.frame("ori" = setdiff(unique(cmc_ori$GENE_NAME), names(ex_by_ge)))
      cmc_genes$new <- cmc_genes$ori
      for (i in 1:nrow(cmc_genes)) {
        incProgress(0, detail = paste0("Update names: ", cmc_genes$ori[i]))
        cmc_genes$new[i] <- UpdateSymbolList(cmc_genes$ori[i])
      }
      cmc_ori$GENE_NAME[cmc_ori$GENE_NAME %in% cmc_genes$ori] <- cmc_genes$new[match(cmc_ori$GENE_NAME[cmc_ori$GENE_NAME %in% cmc_genes$ori], cmc_genes$ori)]
      # remove muts without coordinates
      cmc_ori <- cmc_ori[cmc_ori$`Mutation genome position GRCh37` != "",]  
      # grab seqnames and coordinates
      cmc_ori$chr <- paste0("chr", str_split_fixed(cmc_ori$`Mutation genome position GRCh37`, ":", 2)[,1])
      cmc_ori$chr[cmc_ori$chr == "chr23"] <- "chrX"
      cmc_ori$chr[cmc_ori$chr == "chr24"] <- "chrY"
      cmc_ori$chr[cmc_ori$chr == "chr25"] <- "chrMT"
      cmc_ori$coords <- str_split_fixed(cmc_ori$`Mutation genome position GRCh37`, ":", 2)[,2]
      cmc_ori$start <- as.numeric(str_split_fixed(cmc_ori$coords, "-", 2)[,1])
      cmc_ori$end <- as.numeric(str_split_fixed(cmc_ori$coords, "-", 2)[,2])
      cmc_ori$COSMIC_SAMPLE_POSRATE <- cmc_ori$COSMIC_SAMPLE_MUTATED / cmc_ori$COSMIC_SAMPLE_TESTED
      cmc_ori <- cmc_ori[,c(
        "GENE_NAME","CGC_TIER","Mutation CDS","Mutation AA","chr","start","end","COSMIC_SAMPLE_POSRATE","Mutation Description CDS","Mutation Description AA",
        "GENOMIC_MUTATION_ID","MUTATION_SIGNIFICANCE_TIER"
      )]
      cmc_ori <<- cmc_ori
      incProgress(0, detail = "saving COSMIC db")
      write.fst(cmc_ori, paste0("cmc_", format(Sys.time(), "%Y%m%d_%H%M%S"),".fst"), compress = 100)
      cmc_path <<- list.files(pattern = "cmc_", full.names = T)[length(list.files(pattern = "cmc_"))]
    }
  }
  
  # load and list panel files - need to transition to proper database at some point
  panelfiles <- list.files(path = "panels", pattern = ".panel", full.names = T)
  dbx <- lapply(panelfiles, readRDS)
  dbx_table <- data.frame("panel" = unlist(nameList(dbx, 1)), "sysTimeCreated" = unlist(nameList(dbx, "sysTimeCreated"))) %>%
    group_by(panel) %>%
    mutate("current" = sysTimeCreated == max(sysTimeCreated))
  dbx <- dbx[dbx_table$current]
  
  # Check if PanelCat version is identical across panels, then finalise dbx by naming items, else skip loading and warn
  if (length(unique(unlist(nameList(dbx, "panelCatVersion")))) == 1) {
    names(dbx) <- nameList(dbx, "panelName")
    # how to check if panels with same name are completely different, e.g. compare rows of original target file, TBD
    # get genes of all panels
    all_genes <- vector()
    for (i in 1:length(dbx)) {
      all_genes <- c(all_genes, dbx[[i]][["panelTable"]]$gene)
    }
    all_genes <- data.frame("gene" = sort(unique(all_genes)))
    
    # combine paneltables into one table
    dbx2 <- cbind(data.frame("panel" = rep(dbx[[1]][["panelName"]], nrow(all_genes))),
                  left_join(all_genes, dbx[[1]][["panelTable"]]))
    if (length(dbx) > 1) {
      for (i in 2:length(dbx)) {
        dbx2 <- rbind(dbx2, 
                      cbind(data.frame("panel" = rep(dbx[[i]][["panelName"]], nrow(all_genes))),
                            left_join(all_genes, dbx[[i]][["panelTable"]])))
      }
    }
    status_dbload_success <- T
  } else {
    status_dbload_success <- F
  }
  
  if (length(unique(unlist(nameList(dbx, "txdbv")))) > 1 | 
      length(unique(unlist(nameList(dbx, "clvv")))) > 1 | 
      length(unique(unlist(nameList(dbx, "clvgv")))) > 1 | 
      length(unique(unlist(nameList(dbx, "cmcv")))) > 1) {
    status_dbload_conflict <- T
  } else {
    status_dbload_conflict <- F
  }
  
  # look for existing databases
  txdb_path <- list.files(pattern = "txdb_", full.names = T)[length(list.files(pattern = "txdb_"))]
  cmc_path <- list.files(pattern = "cmc_", full.names = T)[length(list.files(pattern = "cmc_"))]
  clv_path <- list.files(pattern = "clinvar_", full.names = T)[length(list.files(pattern = "clinvar_"))]
  clv_md5_path <- list.files(pattern = "clvmd5", full.names = T)[length(list.files(pattern = "clvmd5"))]
  assRep_path <- list.files(pattern = "assRep_", full.names = T)[length(list.files(pattern = "assRep_"))]
  
  infotext <- readLines("info.txt")
}
# Define UI ---------------------------------------------------------------

css <- "div.dataTables_wrapper  div.dataTables_filter {width: 100%; float: none; text-align: center;}"

ui <- ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(css,"pre { white-space: pre-wrap; word-break: keep-all; }"))),
  headerPanel('PanelCat'),
  
  sidebarPanel(
    conditionalPanel(
      condition = "input.test=='Scatterplot'",
      uiOutput("xpanelsel"),
      uiOutput("xcolsel"),
      uiOutput("ypanelsel"),
      uiOutput("ycolsel"),
      verbatimTextOutput("tableInfo")),
    conditionalPanel(
      condition = "input.test=='Barplot'",
      uiOutput("barpanelsel"),
      selectInput('dtset','Dataset', c("RefSeq","ClinVar","COSMIC")),
      textInput("diff_filter","Diff Filter", value = 0),
      selectInput("vsort", "Sort by", c("gene","minFC","total","ratio","max_cov","max_covp")),
      radioButtons("dtsetrel","display",c("absolute","relative"), selected = "absolute")),
    conditionalPanel(
      condition = "input.test=='Gene'",
      uiOutput("genepanelsel"),
      uiOutput("genesel"),
      selectInput('dtsetgenes','Dataset', c("RefSeq","ClinVar","COSMIC")),
      radioButtons("dtsetrelg","display",c("absolute","relative"), selected = "absolute")),
    conditionalPanel(
      condition = "input.test == 'Exons'",
      uiOutput("expanelsel")),
    conditionalPanel(
      condition = "input.test == 'ExonComp'",
      uiOutput("exoncovpxsel"),
      uiOutput("exoncovpysel"),
      uiOutput("exoncovgsel"),
      uiOutput("exoncovtsel")),
    conditionalPanel(
      condition = "input.test == 'COSMIC'",
      uiOutput("mutpanelsel"),
      checkboxInput("showCmcBl","Hide masked", value = T)),
    conditionalPanel(
      condition = "input.test == 'ClinVar'",
      uiOutput("clvpanelsel"),
      checkboxInput("showClvBl","Hide masked", value = T)),
    conditionalPanel(
      condition = "input.test=='NonCovRate'",
      selectInput("ncovSort","Sort by",c("panel","ncov_all","ncov_all_bl","ncov_targ","ncov_targ_bl"), selected = "panel")),
    conditionalPanel(
      condition = "input.test=='Table'",
      selectInput('tablePanels', 'Panels (choose multiple)', names(dbx), selected = "OFC", multiple = T),
      selectInput('tableVars', 'Variables (choose multiple)', names(dbx2), selected = c("pcb_covp"), multiple = T),
      verbatimTextOutput("tableInfo2")),
    conditionalPanel(# DOES NOT YET PANELS ADDED IN SESSION
      condition = "input.test=='Panels'",
      selectInput('panelInspect', "Select Panel", names(dbx), selected = names(dbx)[1], size = 25, selectize = F)),
    conditionalPanel(
      condition = "input.test=='NewPanel'",
      textInput("pName","Panel Name (abbreviation)"),
      textInput("pfName","Panel Name (full)"),
      fileInput("file1","Panel target file (tab delimited, bed file)"),
      checkboxInput("zeroIndex","Zero-Indexed (.bed convention)?", value = T),
      textInput("pRows","start row, end row"),
      textInput("pCols","Column order: Chr, Start, Stop, Name"),
      fileInput("blackl","Panel blacklist file (optional)"),
      textInput("pmRows","Mask start row, end row"),
      textInput("pmCols","Mask column order: Chr, Start, Stop, Name"),
      actionButton("startb", "Start!"),
      downloadButton("save_state", "Save to file"),
      fileInput("panelUp","processed panel file", accept = ".panel"),
      actionButton("panelUpButton","upload")
      # comment the next line if hosting for others
      ,actionButton("updateb", "UPDATE ALL")
      ),
    width = 3
  ),
  
  mainPanel(
    tabsetPanel(
      id = "test",
      tabPanel("Scatterplot", 
               plotlyOutput("plot1"),
               verbatimTextOutput("stats")),
      tabPanel("Barplot",
               plotOutput("plot2")),
      tabPanel("Gene",
               verbatimTextOutput("unpanels"),
               plotlyOutput("plotGenes")),
      tabPanel("NonCovRate",
               plotlyOutput("nCovPosRate")),
      tabPanel("Table",
               DT::dataTableOutput("table", width = 100)),
      tabPanel("Exons",
               DT::dataTableOutput("exons", width = 100)),
      tabPanel("ExonComp",
               plotOutput("exonCompP")),
      tabPanel("COSMIC",
               DT::dataTableOutput("cmcmuts", width = 100)),
      tabPanel("ClinVar",
               DT::dataTableOutput("clvmuts", width = 100)),
      tabPanel("Panels",
               verbatimTextOutput("viewP")),
      tabPanel("NewPanel", 
               tableOutput("contents"),
               tableOutput("maskFileContent")),
      tabPanel("Info",
               verbatimTextOutput("infotext"))
    ))
)


# # Define MESSY MESSYserver logic  --------------------------------------------------

server <- function(input, output) {
  options(shiny.maxRequestSize=500*1024^2)
  
  # check if DBs loaded
  if (status_dbload_success != T) {
    shinyalert(title = "Error loading panels. Reason: Not all panels analyzed with same Version of PanelCat. Please restore latest functional state, update all, and then add new panels.")
  }
  
  if (status_dbload_conflict == T) {
    shinyalert(title = "Warning! Not all loaded panels were analyzed with identical database versions. Suggest to remove and re-analyze affected panels, or update all.")
  }
  
  output$infotext <- renderText(paste(infotext, collapse = "\n"))
  
  # define reactive input choices
  dbxn <- reactiveValues(panelNames = names(dbx), geneNames = all_genes)
  
  # transition to reactive table later sometime
  #dbxt <- reactiveValues(paneldata = dbx)
  
  # define reactive inputs
  output$xpanelsel <- renderUI({
    selectInput('xpanel', 'X Panel', dbxn$panelNames, selected = dbxn$panelNames[1])
  })
  
  output$ypanelsel <- renderUI({
    selectInput('ypanel', 'Y Panel', dbxn$panelNames, selected = dbxn$panelNames[2])
  })
  
  output$xcolsel <- renderUI({
    selectInput("xcol", "X Variable",names(dbx2)[2:ncol(dbx2)],selected = "pcb_covp")
  })
  
  output$ycolsel <- renderUI({
    selectInput("ycol", "Y Variable",names(dbx2)[2:ncol(dbx2)],selected = "pcb_covp")
  })
  
  output$barpanelsel <- renderUI({
    selectInput('panelbar', 'Panels (choose multiple)', dbxn$panelNames, selected = "OFC", multiple = T)  
  })
  
  output$genepanelsel <- renderUI({
    selectInput('panelgenes', 'Panels (choose multiple)', dbxn$panelNames, selected = dbxn$panelNames, multiple = T)
  })
  
  output$expanelsel <- renderUI({
    selectInput('expanel', 'Panel', dbxn$panelNames, selected = dbxn$panelNames[1])
  })
  
  output$mutpanelsel <- renderUI({
    selectInput('mutpanel', 'Panel', dbxn$panelNames, selected = dbxn$panelNames[2])
  })
  
  output$clvpanelsel <- renderUI({
    selectInput('clvpanel', 'Panel', dbxn$panelNames, selected = "")
  })
  
  output$genesel <- renderUI({
    selectInput('genes', 'Genes (choose multiple)', dbxn$geneNames, selected = "KRAS", multiple = T)
  })
  
  output$exoncovpxsel <- renderUI({
    selectInput('exoncovpx', 'Panel 1', dbxn$panelNames, selected = dbxn$panelNames[1])
  })
  output$exoncovpysel <- renderUI({
    selectInput('exoncovpy', 'Panel 2', dbxn$panelNames, selected = dbxn$panelNames[2])
  })
  output$exoncovgsel <- renderUI({
    selectInput('exoncovg', 'Gene', exoncovd()[["group_name"]], selected = exoncovd()[["group_name"]][1])
  })
  output$exoncovtsel <- renderUI({
    selectInput('exoncovt', 'Transcript', exoncovd2()[["transcript"]], selected = exoncovd()[["transcript"]][1])
  })
  
  output$tableInfo <- renderText(paste("pcb: RefSeq protein coding bases",
                                       "clv: ClinVar variants",
                                       "cmc: CMC mutations",
                                       "tot: total",
                                       "cov: covered (excl. masked)",
                                       "covp: covered % (excl. masked)",
                                       "bl: masked / 'blacklisted'",
                                       "covt: covered (incl. masked)",
                                       "covtp: covered % (incl. masked)",
                                       "ncov: not covered (excl. masked)", sep = "\n \n" ))
  
  output$tableInfo2 <- renderText(paste("pcb: RefSeq protein coding bases",
                                       "clv: ClinVar variants",
                                       "cmc: CMC mutations",
                                       "tot: total",
                                       "cov: covered (excl. masked)",
                                       "covp: covered % (excl. masked)",
                                       "bl: masked / 'blacklisted'",
                                       "covt: covered (incl. masked)",
                                       "covtp: covered % (incl. masked)",
                                       "ncov: not covered (excl. masked)", sep = "\n \n" ))
  
  
  # processed data download -------------------------------------------------
  
  currentpanel <- reactiveVal()
  output$save_state <- downloadHandler(
    filename = function() {
      paste0(currentpanel()[["panelName"]],"_",format(Sys.time(), "%Y%m%d_%H%M%S"), ".panel")
    },
    content = function(file) {
      save_panel <- currentpanel()
      saveRDS(save_panel, file)
    } 
  )  
  
  
  # processed data upload -------------------------------------------------
  
  panelUpFile <- reactive({
    readRDS(input$panelUp$datapath)
  })
  
  observeEvent(input$panelUpButton, {
    if(is.null(input$panelUp$datapath)) {
      shinyalert(title = "No panel file selected for upload!")
    } else {
      if (panelUpFile()[["panelName"]] %in% dbx2$panel) {
        shinyalert(title = "Panel with same name is already loaded")
      } else {
        # update panel dbx
        dbx <<- append(dbx, list(panelUpFile()))
        dbx_table <<- data.frame("panel" = unlist(nameList(dbx, 1)), "sysTimeCreated" = unlist(nameList(dbx, 10))) %>%
          group_by(panel) %>%
          mutate("current" = sysTimeCreated == max(sysTimeCreated))
        dbx <<- dbx[dbx_table$current]
        
        names(dbx) <<- nameList(dbx, "panelName")
        
        # get genes of all panels
        all_genes <<- vector()
        for (i in 1:length(dbx)) {
          all_genes <<- c(all_genes, dbx[[i]][["panelTable"]]$gene)
        }
        all_genes <<- data.frame("gene" = sort(unique(all_genes)))
        
        # combine paneltables into one table
        dbx2 <<- cbind(data.frame("panel" = rep(dbx[[1]][["panelName"]], nrow(all_genes))),
                       left_join(all_genes, dbx[[1]][["panelTable"]]))
        if (length(dbx) > 1) {
          for (i in 2:length(dbx)) {
            dbx2 <<- rbind(dbx2, 
                           cbind(data.frame("panel" = rep(dbx[[i]][["panelName"]], nrow(all_genes))),
                                 left_join(all_genes, dbx[[i]][["panelTable"]])))
          }
        }
        
        dbxn$panelNames <- c(dbxn$panelNames, panelUpFile()[["panelName"]])
        dbxn$geneNames <- all_genes
        shinyalert(title = "complete")
      }
    }
  })
  
  
  
  # X/Y ---------------------------------------------------------------------
  
  # scatter data
  scatter_data <- reactive({
    scatter_dat <- data.frame("x" = dbx2[dbx2$panel == input$xpanel,input$xcol],
                              "y" = dbx2[dbx2$panel == input$ypanel,input$ycol],
                              "gene" = dbx2[dbx2$panel == input$ypanel,"gene"])
    
    scatter_dat$x[(is.na(scatter_dat$x) &  !is.na(scatter_dat$y))] <- 0
    scatter_dat$y[(is.na(scatter_dat$y) &  !is.na(scatter_dat$x))] <- 0
    scatter_dat
  })
  # scatter plots
  output$plot1 <- renderPlotly(
    ggplotly(
      ggplot(d = scatter_data(), aes(text = gene, x = x, y = y)) +
        geom_point(size = 1, alpha = 0.3) +
        theme_minimal() +
        xlab(paste(input$xpanel,input$xcol)) +
        ylab(paste(input$ypanel, input$ycol))
    )
  )
  
  # X/Y overview stats
  setAB <- reactive(
    paste(sort(setdiff(dbx[[input$xpanel]][["panelTable"]][["gene"]],dbx[[input$ypanel]][["panelTable"]][["gene"]])), collapse = " ")
  )
  isectBA <- reactive(
    paste(sort(intersect(dbx[[input$xpanel]][["panelTable"]][["gene"]],dbx[[input$ypanel]][["panelTable"]][["gene"]])), collapse = " ")
  )
  setBA <- reactive(
    paste(sort(setdiff(dbx[[input$ypanel]][["panelTable"]][["gene"]],dbx[[input$xpanel]][["panelTable"]][["gene"]])), collapse = " ")
  )
  
  output$stats <- renderText(paste0(input$xpanel,": ",nrow(as.data.frame(dbx[[input$xpanel]][["panelTable"]]))," genes \n",
                                    input$ypanel,": ",nrow(as.data.frame(dbx[[input$ypanel]][["panelTable"]]))," genes \n\n",
                                    "Genes exclusively in Panel ",input$xpanel,": ",setAB(),"\n\n","Genes exclusively in Panel ",input$ypanel,": ",setBA(),"\n\n","Genes shared: ",isectBA()))
  
  # Panel infos -------------------------------------------------------------
  
  
  output$viewP <- renderText(paste0("Panel Display Name: ", dbx[[input$panelInspect]][["panelName"]],"\n",
                                    "Full Panel Name: ", dbx[[input$panelInspect]][["panelFullName"]],"\n",
                                    "Zero-Based Index: ", dbx[[input$panelInspect]][["zeroIndex"]],"\n",
                                    "Blacklist provided: ", !is.na(dbx[[input$panelInspect]]["blacklist"]),"\n",
                                    "Date analyzed: ", dbx[[input$panelInspect]][["sysTimeCreated"]],"\n",
                                    "RefSeq database: ", dbx[[input$panelInspect]][["txdbv"]],"\n",
                                    "ClinVar database: ", dbx[[input$panelInspect]][["clvv"]],"\n",
                                    "COSMIC CMC database: ", dbx[[input$panelInspect]][["cmcv"]],"\n"))
  
  # panel bars --------------------------------------------------------------
  
  # bar data 
  bars_cov <- reactive({
    bars_cov1 <- dbx2 %>% filter(panel %in% input$panelbar & !is.na(pcb_tot)) %>% filter (case_when(input$dtset == "ClinVar" ~ clv_tot != 0,
                                                                                              input$dtset == "RefSeq" ~ !is.na(pcb_tot),
                                                                                              input$dtset == "COSMIC" ~ cmc_tot != 0)) %>% 
      as_tibble() %>% 
      mutate_at(c("panel", "gene"), as.factor) %>% 
      complete(panel,gene)
    
    bars_cov1 <- left_join(bars_cov1, bars_cov1 %>%
                             replace(is.na(.), 0) %>%
                             #replace(is.nan(.), 0) %>%
                             group_by(gene)%>%
                             summarise(minFC = case_when(input$dtset == "RefSeq" ~ max(pcb_covp) / min(pcb_covp),
                                                         input$dtset == "ClinVar" ~ max(clv_covp) / min(clv_covp),
                                                         input$dtset == "COSMIC" ~ max(cmc_covp) / min(cmc_covp)),
                                       total = case_when(input$dtset == "RefSeq" ~ max(pcb_tot),
                                                         input$dtset == "ClinVar" ~ max(clv_tot),
                                                         input$dtset == "COSMIC" ~ max(cmc_tot)),
                                       ratio = case_when(input$dtset == "RefSeq" ~ pcb_covp[panel == input$panelbar[1]]/pcb_covp[panel == input$panelbar[2]],
                                                         input$dtset == "ClinVar" ~ clv_covp[panel == input$panelbar[1]]/clv_covp[panel == input$panelbar[2]],
                                                         input$dtset == "COSMIC" ~ cmc_covp[panel == input$panelbar[1]]/cmc_covp[panel == input$panelbar[2]]),
                                       max_cov = case_when(input$dtset == "RefSeq" ~ max(pcb_cov),
                                                           input$dtset == "ClinVar" ~ max(clv_cov),
                                                           input$dtset == "COSMIC" ~ max(cmc_cov)),
                                       max_covp = case_when(input$dtset == "RefSeq" ~ max(pcb_covp),
                                                            input$dtset == "ClinVar" ~ max(clv_covp),
                                                            input$dtset == "COSMIC" ~ max(cmc_covp)))) %>%
      filter(minFC > as.numeric(input$diff_filter) | is.nan(minFC))
    
    if (input$vsort != "gene") {
      bars_cov1$gene <- reorder(bars_cov1$gene, bars_cov1[[input$vsort]])
    }
    if (input$vsort == "gene") {
      bars_cov1$gene <- reorder(bars_cov1$gene, rev(order(sort(bars_cov1$gene))))
    }
    
    bars_cov1
  })
  
  bars_tot <- reactive({
    bars_cov() %>%
      replace(is.na(.), 0) %>%
      group_by(gene)%>%
      summarise(width = max(case_when(input$dtset == "RefSeq" ~ pcb_tot,
                                      input$dtset == "ClinVar" ~ clv_tot,
                                      input$dtset == "COSMIC" ~ cmc_tot)))
  })
  
  bars_colvec <- reactive({
    col_vec <- c("grey95",distinctColorPalette(length(input$panelbar)))
    names(col_vec) <- c("total",input$panelbar)
    col_vec
  })
  
  # bar plot
  output$plot2 <- renderPlot({
    validate(
      need(input$panelbar, 'Please select at least one panel!'),
    )
    ggplot() +
        geom_col(d = bars_tot(), aes(x = gene, y = case_when(input$dtsetrel == "absolute" ~ width,
                                                             input$dtsetrel == "relative" ~ 1),
                                     fill = "total")) +
        geom_col(d = bars_cov(), aes(x = gene, y = case_when(input$dtset == "RefSeq" & input$dtsetrel == "absolute" ~ pcb_covt,
                                                             input$dtset == "ClinVar" & input$dtsetrel == "absolute" ~ clv_covt,
                                                             input$dtset == "COSMIC" & input$dtsetrel == "absolute" ~ cmc_covt,
                                                             input$dtset == "RefSeq" & input$dtsetrel == "relative" ~ pcb_covtp,
                                                             input$dtset == "ClinVar" & input$dtsetrel == "relative" ~ clv_covtp,
                                                             input$dtset == "COSMIC" & input$dtsetrel == "relative" ~ cmc_covtp)
                                     , fill = panel), position = "dodge", alpha = 0.3) +
        geom_col(d = bars_cov(), aes(x = gene,  y = case_when(input$dtset == "RefSeq" & input$dtsetrel == "absolute" ~ pcb_cov,
                                                              input$dtset == "ClinVar" & input$dtsetrel == "absolute" ~ clv_cov,
                                                              input$dtset == "COSMIC" & input$dtsetrel == "absolute" ~ cmc_cov,
                                                              input$dtset == "RefSeq" & input$dtsetrel == "relative" ~ pcb_covp,
                                                              input$dtset == "ClinVar" & input$dtsetrel == "relative" ~ clv_covp,
                                                              input$dtset == "COSMIC" & input$dtsetrel == "relative" ~ cmc_covp), 
                                     fill = panel), position = "dodge") +
        ylab(case_when(input$dtset == "RefSeq" & input$dtsetrel == "absolute" ~ "protein coding bases, absolute",
                       input$dtset == "ClinVar" & input$dtsetrel == "absolute" ~ "ClinVar variants, absolute",
                       input$dtset == "COSMIC" & input$dtsetrel == "absolute" ~ "COSMIC mutations, absolute",
                       input$dtset == "RefSeq" & input$dtsetrel == "relative" ~ "protein coding bases, relative",
                       input$dtset == "ClinVar" & input$dtsetrel == "relative" ~ "ClinVar variants, relative",
                       input$dtset == "COSMIC" & input$dtsetrel == "relative" ~ "COSMIC mutations, absolute")) +
        theme_minimal() +
        coord_flip() +
        scale_fill_manual(values = bars_colvec(),
                          guide = guide_legend(reverse=T)) +
        theme(legend.position="right",legend.justification="top")+
      ggtitle("Coverage (opaque) / masked (shaded)")
  }, height = function() {
    10*nrow(bars_cov()) + 200
  })
  
  
  # Genes -------------------------------------------------------------------
  
  # bar data 
  bars_gene <- reactive({
    dbx2[dbx2$panel %in% input$panelgenes & dbx2$gene %in% input$genes & !is.na(dbx2$pcb_tot),] %>% 
      as_tibble() %>% 
      mutate_at(c("panel", "gene"), as.factor) %>% 
      complete(panel,gene)
  })
  
  genes_tot <- reactive({
    bars_gene() %>%
      replace(is.na(.), 0) %>%
      group_by(gene)%>%
      summarise(width = max(case_when(input$dtsetgenes == "RefSeq" & input$dtsetrelg == "absolute" ~ pcb_tot,
                                      input$dtsetgenes == "ClinVar" & input$dtsetrelg == "absolute" ~ clv_tot,
                                      input$dtsetgenes == "COSMIC" & input$dtsetrelg == "absolute" ~ cmc_tot,
                                      input$dtsetrelg == "relative" ~ 1)))
  })
  
  genes_colvec <- reactive({
    col_vec <- c("grey95",distinctColorPalette(length(input$panelgenes)))
    names(col_vec) <- c("total",input$panelgenes)
    col_vec
  })
  
  # gene bar plot
  output$plotGenes <- renderPlotly({
    validate(
      need(input$genes, 'Please select at least one target gene'),
      need(input$panelgenes, 'Please select at least one panel'),
    )
    ggplotly(
      ggplot() +
        geom_col(d = genes_tot(), aes(x = gene, y = width, fill = "total")) +
        geom_col(d = bars_gene(), aes(x = gene, y = case_when(input$dtsetgenes == "RefSeq" & input$dtsetrelg == "absolute" ~ pcb_covt,
                                                              input$dtsetgenes == "ClinVar" & input$dtsetrelg == "absolute" ~ clv_covt,
                                                              input$dtsetgenes == "COSMIC" & input$dtsetrelg == "absolute" ~ cmc_covt,
                                                              input$dtsetgenes == "RefSeq" & input$dtsetrelg == "relative" ~ pcb_covtp,
                                                              input$dtsetgenes == "ClinVar" & input$dtsetrelg == "relative" ~ clv_covtp,
                                                              input$dtsetgenes == "COSMIC" & input$dtsetrelg == "relative" ~ cmc_covtp), 
                                      fill = panel), position = "dodge", alpha = 0.3) +
        geom_col(d = bars_gene(), aes(x = gene,  y = case_when(input$dtsetgenes == "RefSeq" & input$dtsetrelg == "absolute" ~ pcb_cov,
                                                               input$dtsetgenes == "ClinVar" & input$dtsetrelg == "absolute" ~ clv_cov,
                                                               input$dtsetgenes == "COSMIC" & input$dtsetrelg == "absolute" ~ cmc_cov,
                                                               input$dtsetgenes == "RefSeq" & input$dtsetrelg == "relative" ~ pcb_covp,
                                                               input$dtsetgenes == "ClinVar" & input$dtsetrelg == "relative" ~ clv_covp,
                                                               input$dtsetgenes == "COSMIC" & input$dtsetrelg == "relative" ~ cmc_covp), 
                                      fill = panel), position = "dodge") +
        theme_minimal() +
        coord_flip() +
        ggtitle("Coverage (opaque) / masked (shaded)") +
        scale_fill_manual(values = genes_colvec()) %>% reverse_legend_labels(),
      height = case_when((10 * length(input$genes) * length(input$panelgenes)) < 300 ~ 300,
                         (10 * length(input$genes) * length(input$panelgenes)) > 300 ~ (10 * length(input$genes) * length(input$panelgenes)))
    )
  })
  
  unpanelsel <- reactive({
    bars_gene_full <- bars_gene()
    bars_gene_Na <- unique(bars_gene_full$panel[is.na(bars_gene_full$pcb_cov)])
    setdiff(unique(bars_gene_full$panel), bars_gene_Na)
  })
  
  # union panels text
  output$unpanels <-  renderText(c("Panels containing all selected genes: ", paste(unpanelsel(), collapse = " ")))
  
  
  # Non-targeted mutation stats -----------------------------------------------
  
  nCovPosRate_df <- reactive({
    nCovPosRate <- data.table("panel" = unlist(nameList(dbx, "panelName")),
                              "ncov_all_bl" = unlist(nameList(dbx, "cmcNcovPosRateTotal_bl")),
                              "ncov_all" = unlist(nameList(dbx, "cmcNcovPosRateTotal")),
                              "ncov_targ_bl" = unlist(nameList(dbx, "cmcNcovPosRate_bl")), 
                              "ncov_targ" = unlist(nameList(dbx, "cmcNcovPosRate")))
    if (input$ncovSort != "panel") {
      nCovPosRate$panel <- reorder(nCovPosRate$panel, nCovPosRate[[input$ncovSort]])
    }
    if (input$ncovSort == "panel") {
      nCovPosRate$panel <- reorder(nCovPosRate$panel, rev(order(nCovPosRate$panel)))
    }
    melt(nCovPosRate, id.vars = "panel")
  })
  
  output$nCovPosRate <- renderPlotly({
    ggplotly(
      ggplot(d = nCovPosRate_df(), aes(x = panel, y = value, fill = variable)) +
        geom_col(position = position_dodge2(preserve = "single")) +
        theme_minimal() +
        coord_flip()
    )})  
  
  
  # Exon coverage table -----------------------------------------------------

  table_exons <- reactive({
      createExonTable()
    left_join(dbx[[input$expanel]][["exon_coverage"]], ex_by_ge_df1) %>% mutate(covp = cov_width / width)
  })
  
  output$exons <- DT::renderDataTable({
    DT::datatable(table_exons(),
                  filter = list(position = "top", clear = F),
                  options = list(search = list(regex = TRUE, caseInsensitive = T)))
  })
  
#  output$download <- downloadHandler(
#    filename = function(){"thename.csv"}, 
#    content = function(fname){
#      write.csv(thedata(), fname)
#    }
#  )
  

# Exon coverage graph -----------------------------------------------------

  exoncovd <- reactive({
    createExonTable()
  full_join(dbx[[input$exoncovpx]][["exon_coverage"]],
    dbx[[input$exoncovpy]][["exon_coverage"]], 
    by = c("group_name","seqnames","start","end","strand","width"))
  })
  
  exoncovd1 <- reactive({
    exoncovd() %>% filter(group_name == input$exoncovg)
  })
  
  exoncovd2 <- reactive({
    left_join(exoncovd1(), ex_by_ge_df1)
    })
  
  exoncovd3 <- reactive({
    exoncov <- exoncovd2() %>% filter(transcript == input$exoncovt)
    names(exoncov)[names(exoncov) == "cov_width.x"] <- input$exoncovpx
    names(exoncov)[names(exoncov) == "cov_width.y"] <- input$exoncovpy
    exoncov <- melt(exoncov, id.vars = c("group_name","seqnames","start","end","strand","width","transcript","exon_id","exon"))
    exoncov$covp <- exoncov$value / exoncov$width
    exoncov
  })
  
  output$exonCompP <- renderPlot({
    #validate(
    #  need(input$panelbar, 'Please select at least one panel!'),
    #)
    ggplot(d = exoncovd3(), aes(x = variable, y = covp, label = exon)) +
      geom_violin() +
      geom_jitter(position = position_jitter(seed = 1, height = 0)) +
      geom_text_repel(position = position_jitter(seed = 1)) +
      theme_minimal() +
      ylab(paste(input$exoncovg, " exon base coverage, %")) +
      xlab("panel")
  })
  
  # COSMIC coverage table -------------------------------------------------
  
  table_muts <- reactive({
    if (!exists("cmc_ori")) {
      withProgress(message = "Preparation", value = 0, max = 4, {
        loadCosmic(updateDb = FALSE)
      })
    }
    gr_test <- GRanges(dbx[[input$mutpanel]][["panelBed_input"]][["V1"]],
                       IRanges(
                         dbx[[input$mutpanel]][["panelBed_input"]][["V2"]],
                         dbx[[input$mutpanel]][["panelBed_input"]][["V3"]]))
    gr_cmc <- GRanges(cmc_ori$chr, IRanges(cmc_ori$start, cmc_ori$end))
    if (input$showCmcBl == F | length(dbx[[input$mutpanel]][["blacklist"]]) == 1) {
      muts_overlaps <- findOverlaps(gr_test, gr_cmc)
      cmc_ori[muts_overlaps@to]
    } else {
      gr_bl <- GRanges(dbx[[input$mutpanel]][["blacklist"]][["V1"]],
                       IRanges(
                         dbx[[input$mutpanel]][["blacklist"]][["V2"]],
                         dbx[[input$mutpanel]][["blacklist"]][["V3"]]))
      gr_test_bl <- unlist(subtract(gr_test, gr_bl))
      muts_overlaps <- findOverlaps(gr_test_bl, gr_cmc)
      cmc_ori[muts_overlaps@to]
    }
  })
  
  output$cmcmuts <- DT::renderDataTable({
    DT::datatable(table_muts(), 
                  filter = list(position = "top", clear = F),
                  options = list(search = list(regex = TRUE, caseInsensitive = T)))
  })
  
  # CLINVAR coverage table -------------------------------------------------

  
  table_clvmuts <- reactive({
    if (!exists("clinvar")) {
      withProgress(message = "Preparation", value = 0, max = 4, {
        loadClinVar(updateDb = FALSE)
      })
    }
    gr_test <- GRanges(dbx[[input$clvpanel]][["panelBed_input"]][["V1"]],
                       IRanges(
                         dbx[[input$clvpanel]][["panelBed_input"]][["V2"]],
                         dbx[[input$clvpanel]][["panelBed_input"]][["V3"]]))
    gr_clinvar <- GRanges(clinvar$`#CHROM`, IRanges(clinvar$POS, (clinvar$POS+(nchar(clinvar$REF)-1))))
    if (input$showClvBl == F | length(dbx[[input$clvpanel]][["blacklist"]]) == 1) {
      muts_overlaps <- findOverlaps(gr_test, gr_clinvar)
      clinvar[muts_overlaps@to,]
    } else {
      gr_bl <- GRanges(dbx[[input$clvpanel]][["blacklist"]][["V1"]],
                       IRanges(
                         dbx[[input$clvpanel]][["blacklist"]][["V2"]],
                         dbx[[input$clvpanel]][["blacklist"]][["V3"]]))
      gr_test_bl <- unlist(subtract(gr_test, gr_bl))
      muts_overlaps <- findOverlaps(gr_test_bl, gr_clinvar)
      clinvar[muts_overlaps@to,]
    }
  })
  
  output$clvmuts <- DT::renderDataTable({
    DT::datatable(table_clvmuts(),
                  filter = list(position = "top", clear = F),
                  options = list(search = list(regex = TRUE, caseInsensitive = T)))
  })
  
  
  # Panel coverage Tables ------------------------------------------------------------------
  
  table_output <- reactive({
    table <- dcast(as.data.table(dbx2[dbx2$panel %in% input$tablePanels,]), gene~panel, value.var = input$tableVars)
    table <- table[apply(table, 1, count_nna_func) > 1,]
    table
  })
  
  output$table <- DT::renderDataTable({
    validate(
      need(input$tablePanels, 'Select at least one panel'),
      need(input$tableVars, "Select at least one variable")
    )
    DT::datatable(table_output())
  })
  
  
  # NEW PANEL  ---------------------------------------------------------
  
  # target file
  panelInput <- reactiveValues(data=NULL, mask=NULL)
  
  observe({
    req(input$file1)
    panelInput$data <- read.table(input$file1$datapath, fill = T, row.names = NULL)
  })
  
  rowVec <- reactive({
    rowVecs <- as.vector(as.numeric(str_split_fixed(input$pRows,",",2)))
    if (is.na(rowVecs[1])) {
      rowVecs[1] <- 1
    }
    if (is.na(rowVecs[2])) {
      rowVecs[2] <- nrow(panelInput$data)
    }
    seq(rowVecs[1],rowVecs[2])
  })
  
  colVec <- reactive({
    colVecs <- as.vector(as.numeric(str_split_fixed(input$pCols,",",3)))
    if (sum(is.na(colVecs)) > 0) {
      colVecs <- seq(1,ncol(panelInput$data))
    }
    if (sum(is.na(colVecs)) == 0) {
      colVecs <- colVecs
    }
    colVecs
  })
  
  output$contents <- renderTable({
    validate(
      need(!is.null(panelInput$data), 'Please provide a target region file!'),
    )
    rbind(head(panelInput$data[rowVec(),colVec()], 3), tail(panelInput$data[rowVec(),colVec()],3))
  })
  
  # Mask file
  observe({
    req(input$blackl)
    panelInput$mask <- read.table(input$blackl$datapath, fill = T, row.names = NULL)
  })
  
  mrowVec <- reactive({
    rowVecs <- as.vector(as.numeric(str_split_fixed(input$pmRows,",",2)))
    if (is.na(rowVecs[1])) {
      rowVecs[1] <- 1
    }
    if (is.na(rowVecs[2])) {
      rowVecs[2] <- nrow(panelInput$mask)
    }
    seq(rowVecs[1],rowVecs[2])
  })
  
  mcolVec <- reactive({
    colVecs <- as.vector(as.numeric(str_split_fixed(input$pmCols,",",3)))
    if (sum(is.na(colVecs)) > 0) {
      colVecs <- seq(1,ncol(panelInput$mask))
    }
    if (sum(is.na(colVecs)) == 0) {
      colVecs <- colVecs
    }
    colVecs
  })
  
  output$maskFileContent <- renderTable({
    validate(
      need(!is.null(panelInput$mask), "Please provide a mask file (optional!)"),
    )
    rbind(head(panelInput$mask[mrowVec(),mcolVec()], 3), tail(panelInput$mask[mrowVec(),mcolVec()],3))
  })
  
  observeEvent(input$startb, {
    
    #reorder file columns
    trget <- panelInput$data[rowVec(),colVec()]
    names(trget) <- c("V1", "V2","V3")
    
    # apply zero-based index correction if dealing with .bed convention
    if (input$zeroIndex == T) {
      trget$V2 <- trget$V2 + 1
    }
    
    if (is.null(panelInput$mask)) {
      blacklist <- NA
      gr_blacklist <- GRanges()
      blacklist_check <- T
    } else if (!is.null(panelInput$mask)) {
      blacklist <- panelInput$mask[mrowVec(),mcolVec()]
      names(blacklist) <- c("V1", "V2","V3")
      if ((sum(anyNA(as.numeric(blacklist$V2)), anyNA(as.numeric(blacklist$V3))) != 0)) {
        blacklist_check = F
      } else {
        blacklist_check <- T
      }
    }
    
    
    if (input$pName %in% dbx2$panel) {
      shinyalert(title = "Panel Name already in use. Please provide a new unique name.")
    } else if (is.null(input$pName)) {
      shinyalert(title = "Please provide a new unique panel name.")
    } else if ((sum(anyNA(as.numeric(trget$V2)), anyNA(as.numeric(trget$V3))) != 0)) {
      shinyalert(title = "Specified panel file columns contain non-numeric entries. Please review.")
    } else if (blacklist_check == F) {
      shinyalert(title = "Specified mask file columns contain non-numeric entries. Please review.")
    } else {

      withProgress(message = "", value = 0, max = 8, {
        # load databases
        loadRefSeq(updateDb = FALSE)
        loadClinVar(updateDb = FALSE)
        loadCosmic(updateDb = FALSE)
        cmc <- cmc_ori[cmc_ori$MUTATION_SIGNIFICANCE_TIER != "Other",]
        
        # make Granges of panel
        trget$V2 <- as.numeric(trget$V2)
        trget$V3 <- as.numeric(trget$V3)
        gr_test <- GRanges(trget$V1, IRanges(trget$V2, trget$V3))
        
        # make blacklist ranges
        if (!is.null(panelInput$mask)) {
        blacklist$V2 <- as.numeric(blacklist$V2)
        blacklist$V3 <- as.numeric(blacklist$V3)
        names(blacklist) <- c("V1", "V2","V3")
        gr_blacklist <- GRanges(blacklist$V1, IRanges(blacklist$V2, blacklist$V3))
        }
        
        gr_test_bl <- reduce(unlist(subtract(gr_test, gr_blacklist)))
        
        #find overlapping exons in txdb
        incProgress(1, detail = "Find targeted exons")
        ex_by_ge_overlap <- findOverlaps(gr_test, unlist(ex_by_ge)) 
        ex_by_ge_overlap <- unlist(ex_by_ge)[ex_by_ge_overlap@to]
        gl_exp <<- data.frame("refseq" = unique(names(ex_by_ge_overlap)))
        
        # identify unique exons and their coverage
        
        incProgress(1, detail = "Find unique exon coverage (optional)")
        gr_test_red <- reduce(gr_test)
        exons <- ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq]
        ex_fo <- findOverlaps(gr_test_red, exons)
        ex_by_ge_all_df <- subset(
          left_join(as.data.table(exons),
                    as.data.table(pintersect(gr_test_red[queryHits(ex_fo)], exons[subjectHits(ex_fo)])) %>% 
                      filter(hit == T) %>%
                      group_by(exon_id) %>%
                      summarise(cov_width = sum(width)),
                    by = "exon_id"),
          select = c(group_name, seqnames, start, end, strand, width, cov_width)) %>%
          distinct() %>%
          replace(is.na(.), 0)
        
        # reduce exon ranges
        incProgress(1, detail = "Reducing exon GRanges")
        exons_red <- reduce(ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq])
        
        ex_fo1 <- findOverlaps(gr_test_red, exons_red)
        ex_fo_pint1 <- pintersect(gr_test_red[queryHits(ex_fo1)], exons_red[subjectHits(ex_fo1)])
        table_refseq_cov <- as.data.table(ex_fo_pint1) %>% group_by(group_name) %>% summarise(pcb_covt = sum(width))
        
        ex_fo2 <- findOverlaps(gr_test_bl, exons_red)
        ex_fo_pint2 <- pintersect(gr_test_bl[queryHits(ex_fo2)], exons_red[subjectHits(ex_fo2)])
        table_refseq_covbl <- as.data.table(ex_fo_pint2) %>% group_by(group_name) %>% summarise(pcb_cov = sum(width))
        
        ex_fo3 <- findOverlaps(gr_blacklist, exons_red)
        ex_fo_pint3 <- pintersect(gr_blacklist[queryHits(ex_fo3)], exons_red[subjectHits(ex_fo3)])
        table_refseq_bl <- as.data.table(ex_fo_pint3) %>% group_by(group_name) %>% summarise(pcb_bl = sum(width))
        
        table_refseq <<- left_join(
          left_join(table_refseq_cov, table_refseq_bl),
          table_refseq_covbl) %>% mutate(pcb_tot = sum(width(exons_red)), pcb_ncov = pcb_tot - pcb_covt, pcb_covp = pcb_cov / pcb_tot, pcb_covtp = pcb_covt / pcb_tot)
        colnames(table_refseq)[1] <- "gene"
        table_refseq <- table_refseq[,c("gene","pcb_ncov","pcb_covt","pcb_bl","pcb_cov","pcb_tot","pcb_covp","pcb_covtp")]
        
        # CLINVAR TABLE
        incProgress(1, detail = "Create ClinVar table")
        
        # filter clinvar to genes of test
        clinvar_test <- clinvar[clinvar$gene %in% gl_exp$refseq,]
        
        # filter clinvar to pathogenic and likely pathogenic
        clinvar_select <- clinvar_test[str_starts(clinvar_test$clnsig, "Pathogenic") | str_starts(clinvar_test$clnsig, "Likely_pathogenic"),]
        table_clinvar_tot <- as.data.frame(table(clinvar_select$gene))
        
        # create clinvar ranges and overlap with test ranges
        gr_clinvar <- GRanges(clinvar_select$`#CHROM`, IRanges(clinvar_select$POS, (clinvar_select$POS+(nchar(clinvar_select$REF)-1))), 
                               REF = clinvar_select$REF, ALT = clinvar_select$ALT)
        
        # covered variants (excluding blacklisted)
        clinvar_fo_bl <- findOverlaps(gr_clinvar, gr_test_bl, type = "within")
        table_clinvar_covbl <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo_bl),])[["gene"]]))
        
        # total and blacklisted variants
        clinvar_fo <- findOverlaps(gr_clinvar, gr_test_red, type = "within")
        table_clinvar_covt <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo),])[["gene"]]))
        
        clinvar_fo_blonly <- findOverlaps(gr_clinvar, gr_blacklist, type = "within")
        if(!is.null(panelInput$mask)) {
          table_clinvar_bl <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo_blonly),])[["gene"]]))
        } else {
          table_clinvar_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }
        
        # NON covered variants (excluding blacklisted)
        #clinvar_ncov <- distinct(clinvar_select[-unique(queryHits(clinvar_fo)),])[["gene"]]
        table_clinvar_ncov <- as.data.frame(table(distinct(clinvar_select[-unique(queryHits(clinvar_fo)),])[["gene"]]))
        
        # join tables
        table_clinvar <<- left_join(
          left_join(
            left_join(
              left_join(
                full_join(data.frame("Var1" = gl_exp$refseq), table_clinvar_ncov, by = "Var1"),
                table_clinvar_covt, by = "Var1"),
              table_clinvar_bl, by = "Var1"),
            table_clinvar_covbl, by = "Var1"),
          table_clinvar_tot, by = "Var1")
        colnames(table_clinvar) <- c("gene","clv_ncov","clv_covt","clv_bl","clv_cov","clv_tot")
        table_clinvar$clv_covp <- table_clinvar$clv_cov/table_clinvar$clv_tot
        table_clinvar$clv_covtp <- table_clinvar$clv_covt/table_clinvar$clv_tot
        table_clinvar[is.na(table_clinvar)] <- 0
        
        # Cosmic ------------------------------------------------------------------
        incProgress(1, detail = "Create COSMIC table")
        # make ranges, find overlap (blacklisted)
        gr_cmc <<- GRanges(cmc$chr, IRanges(cmc$start, cmc$end))
        
        # total muts in targeted genes
        cmc_test <<- cmc[cmc$GENE_NAME %in% gl_exp$refseq]
        table_cmc_tot <<- as.data.frame(table(cmc_test$GENE_NAME))
        
        # covered variants (excluding blacklisted)
        cmc_fo_bl <<- findOverlaps(gr_cmc, gr_test_bl, type = "within")
        cmc_cov_bl <<- distinct(cmc[cmc_fo_bl@from,])
        table_cmc_covbl <<- as.data.frame(sort(table(cmc_cov_bl$GENE_NAME)))
        
        # total and explicitly blacklisted variants
        cmc_fo <<- findOverlaps(gr_cmc, gr_test_red, type = "within")
        cmc_cov <<- distinct(cmc[cmc_fo@from,])
        table_cmc_covt <<- as.data.frame(table(cmc_cov$GENE_NAME))
        cmc_bl <<- cmc_cov[!(cmc_cov$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        table_cmc_bl <<- as.data.frame(table(cmc_bl$GENE_NAME))
        if(!is.null(panelInput$mask)) {
          table_cmc_bl <<- as.data.frame(table(cmc_bl$GENE_NAME))
        } else {
          table_cmc_bl <<- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }
        
        # NON covered variants (excluding blacklisted)
        cmc_ncov <<- cmc_test[!(cmc_test$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
        if (nrow(cmc_ncov) > 0) {
          table_cmc_ncov <<- as.data.frame(table(cmc_ncov$GENE_NAME))
        } else {
          table_cmc_ncov <<- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }
        
        # NON covered including blacklisted
        cmc_ncovbl <<- cmc_test[!(cmc_test$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        #table_cmc_ncovbl <<- as.data.frame(table(cmc_ncovbl$GENE_NAME))
        
        # rate of non-covered mutations - only in panel target genes
        cmc_ncov_posTestRate <<- sum(cmc_ncov$COSMIC_SAMPLE_POSRATE)
        cmc_ncovbl_posTestRate <<- sum(cmc_ncovbl$COSMIC_SAMPLE_POSRATE)
        
        # all non-covered mutations (including non-target genes)
        cmc_ncov_total <<- cmc[!(cmc$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
        cmc_ncovbl_total <<- cmc[!(cmc$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        
        # rate of non-covered mutations in all genes (including non targeted genes)
        cmc_ncov_posTestRateTotal <<- sum(cmc_ncov_total$COSMIC_SAMPLE_POSRATE)
        cmc_ncovbl_posTestRateTotal <<- sum(cmc_ncovbl_total$COSMIC_SAMPLE_POSRATE)
        
        # join tables
        table_cmc <<- left_join(
          left_join(
            left_join(
              left_join(
                full_join(data.frame("Var1" = gl_exp$refseq), table_cmc_ncov, by = "Var1"),
                table_cmc_covt, by = "Var1"),
              table_cmc_bl, by = "Var1"),
            table_cmc_covbl, by = "Var1"),
          table_cmc_tot, by = "Var1")
        colnames(table_cmc) <- c("gene","cmc_ncov","cmc_covt","cmc_bl","cmc_cov","cmc_tot")
        table_cmc$cmc_covp <- table_cmc$cmc_cov/table_cmc$cmc_tot
        table_cmc$cmc_covtp <- table_cmc$cmc_covt/table_cmc$cmc_tot
        table_cmc[is.na(table_cmc)] <- 0
        
        # list panel items
        panel <- list(panelName = input$pName, 
                      panelFullName = input$pfName,
                      panelTable = as_tibble(full_join(full_join(table_refseq, table_cmc, by = "gene"), table_clinvar, by = "gene")),
                      panelBed_ori = panelInput$data,
                      panelBed_input = trget,
                      zeroIndex = input$zeroIndex,
                      blacklist = blacklist,
                      cmcNcovPosRate = cmc_ncov_posTestRate,
                      cmcNcovPosRate_bl = cmc_ncovbl_posTestRate, 
                      cmcNcovPosRateTotal = cmc_ncov_posTestRateTotal,
                      cmcNcovPosRateTotal_bl = cmc_ncovbl_posTestRateTotal,
                      sysTimeCreated = Sys.time(),
                      exon_coverage = ex_by_ge_all_df,
                      txdbv = txdb_path,
                      clvv = clv_path,
                      cmcv = cmc_path,
                      panelCatVersion = PanCatv
        )
        # comment the next line if hosting for others
        saveRDS(panel, paste0("panels/", panel[["panelName"]],"_",format(Sys.time(), "%Y%m%d_%H%M%S"), ".panel"))
        
        currentpanel(panel)
        
        # update panel dbx
        dbx <<- append(dbx, list(panel))
        dbx_table <<- data.frame("panel" = unlist(nameList(dbx, 1)), "sysTimeCreated" = unlist(nameList(dbx, 10))) %>%
          group_by(panel) %>%
          mutate("current" = sysTimeCreated == max(sysTimeCreated))
        dbx <<- dbx[dbx_table$current]
        
        names(dbx) <<- nameList(dbx, "panelName")
        
        # how to check if panels with same name are completely different
        # compare rows of original target file
        # TBD
        
        # get genes of all panels
        all_genes <<- vector()
        for (i in 1:length(dbx)) {
          all_genes <<- c(all_genes, dbx[[i]][["panelTable"]]$gene)
        }
        all_genes <<- data.frame("gene" = sort(unique(all_genes)))
        
        # combine paneltables into one table
        dbx2 <<- cbind(data.frame("panel" = rep(dbx[[1]][["panelName"]], nrow(all_genes))),
                       left_join(all_genes, dbx[[1]][["panelTable"]]))
        if (length(dbx) > 1) {
          for (i in 2:length(dbx)) {
            dbx2 <<- rbind(dbx2, 
                           cbind(data.frame("panel" = rep(dbx[[i]][["panelName"]], nrow(all_genes))),
                                 left_join(all_genes, dbx[[i]][["panelTable"]])))
          }
        }
        
        dbxn$panelNames <- c(dbxn$panelNames, input$pName)
        dbxn$geneNames <- all_genes
        reset("pName")
        reset("pfName")
        reset("file1")
        reset("pRows")
        reset("pCols")
        reset("blackl")
        panelInput$data <- NULL
        panelInput$mask <- NULL
        shinyalert(title = "complete")
      })
    }
  })
  
  
  # UPDATE ALL --------------------------------------------------------------
  
  observeEvent(input$updateb, {
    withProgress(message = "Preparation", value = 0, max = 4, {
      
      loadRefSeq(updateDb = TRUE)
      loadClinVar(updateDb = updateCheckClv())
      loadCosmic(updateDb = FALSE)
      cmc <- cmc_ori[cmc_ori$MUTATION_SIGNIFICANCE_TIER != "Other",]
    })
    
    withProgress(message = "Updating", value = 0, max = length(dbx), {
      for (j in 1:length(dbx)) {
        incProgress(1, detail = dbx[[j]][["panelName"]])
        
        
        # PROCESS BED FILE
        #reorder file columns
        trget <- dbx[[j]][["panelBed_input"]]
        
        # make Granges of panel
        gr_test <- GRanges(trget$V1, IRanges(trget$V2, trget$V3))
        
        # blacklist
        blacklist <- dbx[[j]][["blacklist"]]
        
        if (length(dbx[[j]][["blacklist"]]) == 1) {
          gr_blacklist <- GRanges()
        }
        if (length(dbx[[j]][["blacklist"]]) >= 3) {
          gr_blacklist <- GRanges(blacklist$V1, IRanges(blacklist$V2, blacklist$V3))
        }
        gr_test_bl <- reduce(unlist(subtract(gr_test, gr_blacklist)))
        
        #find overlapping exons in txdb
        ex_by_ge_overlap <- findOverlaps(gr_test, unlist(ex_by_ge)) 
        ex_by_ge_overlap <- unlist(ex_by_ge)[ex_by_ge_overlap@to]
        gl_exp <- data.frame("refseq" = unique(names(ex_by_ge_overlap)))
        
        # identify unique exons and their coverage
        gr_test_red <- reduce(gr_test)
        exons <- ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq]
        ex_fo <- findOverlaps(gr_test_red, exons)
        ex_by_ge_all_df <- subset(
          left_join(as.data.table(exons),
                    as.data.table(pintersect(gr_test_red[queryHits(ex_fo)], exons[subjectHits(ex_fo)])) %>% 
                      filter(hit == T) %>%
                      group_by(exon_id) %>%
                      summarise(cov_width = sum(width)),
                    by = "exon_id"),
          select = c(group_name, seqnames, start, end, strand, width, cov_width)) %>%
          distinct() %>%
          replace(is.na(.), 0)
        
        # reduce exon ranges
        exons_red <- reduce(ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq])
        
        ex_fo1 <- findOverlaps(gr_test_red, exons_red)
        ex_fo_pint1 <- pintersect(gr_test_red[queryHits(ex_fo1)], exons_red[subjectHits(ex_fo1)])
        table_refseq_cov <- as.data.table(ex_fo_pint1) %>% group_by(group_name) %>% summarise(pcb_covt = sum(width))
        
        ex_fo2 <- findOverlaps(gr_test_bl, exons_red)
        ex_fo_pint2 <- pintersect(gr_test_bl[queryHits(ex_fo2)], exons_red[subjectHits(ex_fo2)])
        table_refseq_covbl <- as.data.table(ex_fo_pint2) %>% group_by(group_name) %>% summarise(pcb_cov = sum(width))
        
        ex_fo3 <- findOverlaps(gr_blacklist, exons_red)
        ex_fo_pint3 <- pintersect(gr_blacklist[queryHits(ex_fo3)], exons_red[subjectHits(ex_fo3)])
        table_refseq_bl <- as.data.table(ex_fo_pint3) %>% group_by(group_name) %>% summarise(pcb_bl = sum(width))
        
        table_refseq <- left_join(
          left_join(table_refseq_cov, table_refseq_bl),
          table_refseq_covbl) %>% mutate(pcb_tot = sum(width(exons_red)), pcb_ncov = pcb_tot - pcb_covt, pcb_covp = pcb_cov / pcb_tot, pcb_covtp = pcb_covt / pcb_tot)
        colnames(table_refseq)[1] <- "gene"
        table_refseq <- table_refseq[,c("gene","pcb_ncov","pcb_covt","pcb_bl","pcb_cov","pcb_tot","pcb_covp","pcb_covtp")]
        
        # CLINVAR TABLE

        # filter clinvar to genes of test
        clinvar_test <- clinvar[clinvar$gene %in% gl_exp$refseq,]
        
        # filter clinvar to pathogenic and likely pathogenic
        clinvar_select <- clinvar_test[str_starts(clinvar_test$clnsig, "Pathogenic") | str_starts(clinvar_test$clnsig, "Likely_pathogenic"),]
        table_clinvar_tot <- as.data.frame(table(clinvar_select$gene))
        
        # create clinvar ranges and overlap with test ranges
        gr_clinvar <- GRanges(clinvar_select$`#CHROM`, IRanges(clinvar_select$POS, (clinvar_select$POS+(nchar(clinvar_select$REF)-1))), 
                               REF = clinvar_select$REF, ALT = clinvar_select$ALT)
        
        # covered variants (excluding blacklisted)
        clinvar_fo_bl <- findOverlaps(gr_clinvar, gr_test_bl, type = "within")
        table_clinvar_covbl <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo_bl),])[["gene"]]))
        
        # total and blacklisted variants
        clinvar_fo <- findOverlaps(gr_clinvar, gr_test_red, type = "within")
        table_clinvar_covt <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo),])[["gene"]]))
        
        clinvar_fo_blonly <- findOverlaps(gr_clinvar, gr_blacklist, type = "within")
        if(sum(is.na(blacklist)) == 0) {
          table_clinvar_bl <- as.data.frame(table(distinct(clinvar_select[queryHits(clinvar_fo_blonly),])[["gene"]]))
        } else {
          table_clinvar_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }
        
        # NON covered variants (excluding blacklisted)
        #clinvar_ncov <- distinct(clinvar_select[-unique(queryHits(clinvar_fo)),])[["gene"]]
        table_clinvar_ncov <- as.data.frame(table(distinct(clinvar_select[-unique(queryHits(clinvar_fo)),])[["gene"]]))
        
        # join tables
        table_clinvar <- left_join(
          left_join(
            left_join(
              left_join(
                full_join(data.frame("Var1" = gl_exp$refseq), table_clinvar_ncov, by = "Var1"),
                table_clinvar_covt, by = "Var1"),
              table_clinvar_bl, by = "Var1"),
            table_clinvar_covbl, by = "Var1"),
          table_clinvar_tot, by = "Var1")
        colnames(table_clinvar) <- c("gene","clv_ncov","clv_covt","clv_bl","clv_cov","clv_tot")
        table_clinvar$clv_covp <- table_clinvar$clv_cov/table_clinvar$clv_tot
        table_clinvar$clv_covtp <- table_clinvar$clv_covt/table_clinvar$clv_tot
        table_clinvar[is.na(table_clinvar)] <- 0
        
        # Cosmic ------------------------------------------------------------------
        # make ranges, find overlap (blacklisted)
        gr_cmc <- GRanges(cmc$chr, IRanges(cmc$start, cmc$end))
        
        # total muts in targeted genes
        cmc_test <- cmc[cmc$GENE_NAME %in% gl_exp$refseq]
        table_cmc_tot <- as.data.frame(table(cmc_test$GENE_NAME))
        
        # covered variants (excluding blacklisted)
        cmc_fo_bl <- findOverlaps(gr_cmc, gr_test_bl, type = "within")
        cmc_cov_bl <- distinct(cmc[cmc_fo_bl@from,])
        table_cmc_covbl <- as.data.frame(sort(table(cmc_cov_bl$GENE_NAME)))
        
        # total and explicitly blacklisted variants
        cmc_fo <- findOverlaps(gr_cmc, gr_test_red, type = "within")
        cmc_cov <- distinct(cmc[cmc_fo@from,])
        table_cmc_covt <- as.data.frame(table(cmc_cov$GENE_NAME))
        cmc_bl <- cmc_cov[!(cmc_cov$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
        if(sum(is.na(blacklist)) == 0) {
          table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
        } else {
          table_cmc_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }
        
        # NON covered variants (excluding blacklisted)
        cmc_ncov <- cmc_test[!(cmc_test$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
        if (nrow(cmc_ncov) > 0) {
          table_cmc_ncov <<- as.data.frame(table(cmc_ncov$GENE_NAME))
        } else {
          table_cmc_ncov <<- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
        }

        # NON covered including blacklisted
        cmc_ncovbl <- cmc_test[!(cmc_test$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        #table_cmc_ncovbl <- as.data.frame(table(cmc_ncovbl$GENE_NAME))
        
        # rate of non-covered mutations - only in panel target genes
        cmc_ncov_posTestRate <- sum(cmc_ncov$COSMIC_SAMPLE_POSRATE)
        cmc_ncovbl_posTestRate <- sum(cmc_ncovbl$COSMIC_SAMPLE_POSRATE)
        
        # all non-covered mutations (including non-target genes)
        cmc_ncov_total <- cmc[!(cmc$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
        cmc_ncovbl_total <- cmc[!(cmc$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
        
        # rate of non-covered mutations in all genes (including non targeted genes)
        cmc_ncov_posTestRateTotal <- sum(cmc_ncov_total$COSMIC_SAMPLE_POSRATE)
        cmc_ncovbl_posTestRateTotal <- sum(cmc_ncovbl_total$COSMIC_SAMPLE_POSRATE)
        
        # join tables
        table_cmc <- left_join(
          left_join(
            left_join(
              left_join(
                full_join(data.frame("Var1" = gl_exp$refseq), table_cmc_ncov, by = "Var1"),
                table_cmc_covt, by = "Var1"),
              table_cmc_bl, by = "Var1"),
            table_cmc_covbl, by = "Var1"),
          table_cmc_tot, by = "Var1")
        colnames(table_cmc) <- c("gene","cmc_ncov","cmc_covt","cmc_bl","cmc_cov","cmc_tot")
        table_cmc$cmc_covp <- table_cmc$cmc_cov/table_cmc$cmc_tot
        table_cmc$cmc_covtp <- table_cmc$cmc_covt/table_cmc$cmc_tot
        table_cmc[is.na(table_cmc)] <- 0
        
        # list panel items
        panel <- list(panelName = dbx[[j]][["panelName"]], 
                      panelFullName = dbx[[j]][["panelFullName"]],
                      panelTable = as_tibble(full_join(full_join(table_refseq, table_cmc, by = "gene"), table_clinvar, by = "gene")),
                      panelBed_ori = dbx[[j]][["panelBed_ori"]],
                      panelBed_input = trget,
                      zeroIndex = dbx[[j]][["zeroIndex"]],
                      blacklist = blacklist,
                      cmcNcovPosRate = cmc_ncov_posTestRate,
                      cmcNcovPosRate_bl = cmc_ncovbl_posTestRate, 
                      cmcNcovPosRateTotal = cmc_ncov_posTestRateTotal,
                      cmcNcovPosRateTotal_bl = cmc_ncovbl_posTestRateTotal,
                      sysTimeCreated = Sys.time(),
                      exon_coverage = ex_by_ge_all_df,
                      txdbv = txdb_path,
                      clvv = clv_path,
                      cmcv = cmc_path,
                      panelCatVersion = PanCatv
        )
        saveRDS(panel, paste0("panels/", panel[["panelName"]],"_",format(Sys.time(), "%Y%m%d_%H%M%S"), ".panel"))
      }
    })
    shinyalert(title = "complete! Please restart application.")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
