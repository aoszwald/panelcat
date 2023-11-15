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
  if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)}
  if(!require(randomcoloR)){
    install.packages("randomcoloR")
    library(randomcoloR)}
  if(!require(RSQLite)){
    install.packages("RSQLite")
    library(RSQLite)}
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")} 
  if(!require(GenomicFeatures)){
    BiocManager::install("GenomicFeatures")
    library(GenomicFeatures)}
  if(!require(Seurat)){
    BiocManager::install("Seurat")}

  options(repos = BiocManager::repositories())
  options(timeout = 360)
  
  # Functions ---------------------------------------------------------------
  {
    # set true to always skip database update
    debugMode <- F
    
    # PanCatversion
    PanCatv <- 29
    
    count_nna_func <- function(x) sum(!is.na(x))
    
    # reverse plotly legend labels
    reverse_legend_labels <- function(plotly_plot) {
      n_labels <- length(plotly_plot$x$data)
      plotly_plot$x$data[1:n_labels] <- plotly_plot$x$data[n_labels:1]
      plotly_plot
    }
  
    remove_alpha_legend_genes_reverse <- function(plotly_plot) {
      plotly_plot$x$data[[length(plotly_plot$x$data)]]$showlegend <- T
      panelsChecked <- plotly_plot$x$data[[length(plotly_plot$x$data)]]$name
      for (i in (length(plotly_plot$x$data)-1):1) {
        if (!(plotly_plot$x$data[[i]]$name %in% panelsChecked)) {
          plotly_plot$x$data[[i]]$showlegend <- T
          panelsChecked <- c(panelsChecked, plotly_plot$x$data[[i]]$name)
        } else {
          plotly_plot$x$data[[i]]$showlegend <- F
        }
      }
      panelOrder <- plotly_plot$x$data[[1]]$name
      for (i in 2: length(plotly_plot$x$data)) {
        panelOrder <- c(panelOrder, plotly_plot$x$data[[i]]$name)
      }
      panelOrderu <- panelOrder[!duplicated(panelOrder)]
      plotly_plot$x$data[2:length(plotly_plot$x$data)] <- plotly_plot$x$data[rev(order(match(panelOrder, panelOrderu))[2:length(panelOrder)])]
      plotly_plot
    }
    
    remove_alpha_legend <- function(plotly_plot, inputP) {
      for (i in 1:length(plotly_plot$x$data)) {
        if (between(i, 2, (length(inputP)+1))) {
          plotly_plot$x$data[[i]]$showlegend <- F
        } else {
          plotly_plot$x$data[[i]]$showlegend <- T
        }
      }
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
    prepRefSeq <- function(updateDb) {
      withProgress(message = "Preparing RefSeq database", {
        # if database path does not exist, create new
        if ((length(txdb_path) == 0 | updateDb == T) & debugMode == F) {
          incProgress(0.3, detail = "Downloading")
          txdb <- makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz",
                                   dataSource = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz")
          incProgress(0.3, detail = "Processing")
          assRep <- read.table("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt", fill = T)
          seqlevels(txdb) <- seqlevels(txdb)[1:24]
          seqlevels(txdb) <- assRep$V11[1:24]
          ex_by_ge <- exonsBy(txdb, by="gene")
          incProgress(0.3, detail = "Saving")
          saveRDS(ex_by_ge, paste0("txdb_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
        }
        # if exon table does not exist, or if Db was updated, create
        if ((length(txdb_path) == 0 | length(exdb_path) == 0 | updateDb == T)& debugMode == F) {
          incProgress(0.3, detail = "Create exon table")
          ex_by_ge_df1 <- (as.data.table(unlist(ex_by_ge)))
          splitnames <- str_split_fixed(ex_by_ge_df1$exon_name, "-",3)
          colnames(splitnames) <- c("name","transcript","exon")
          ex_by_ge_df1 <- cbind(ex_by_ge_df1, splitnames)
          ex_by_ge_df1 <- subset(ex_by_ge_df1, select = -c(exon_name,name))
          ex_by_ge_df1$exon <- as.numeric(ex_by_ge_df1$exon)
          # save exon table and path for later
          saveRDS(ex_by_ge_df1, paste0("exdb_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
          txdb_path <<- list.files(pattern = "txdb_", full.names = T)[length(list.files(pattern = "txdb_"))]
          exdb_path <<- list.files(pattern = "exdb_", full.names = T)[length(list.files(pattern = "exdb_"))]
        }
      })
    }
    
    loadRefSeq <- function(force) {
      prepRefSeq(updateDb = F)
      if (exists("ex_by_ge") == F | force == T) {
        withProgress(message = "Loading RefSeq database", {
          ex_by_ge <<- readRDS(txdb_path)
        })
      } 
    }
    
    loadExDb <- function() {
      if (length(exdb_path) == 0) 
        prepRefSeq(updateDb = F)
      if (exists("ex_by_ge_df1") == F) {
        ex_by_ge_df1 <<- readRDS(exdb_path)
      }
    }
    
    # CLINVAR update check (only used via "update all")
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
    
    # Clinvar loader / creator
    prepClinvar <- function(updateDb) {
      withProgress(message = "Preparing ClinVar database", {
        # first, check if db exists or forced update, if not, create
        if ((length(clv_path) == 0 | updateDb == T) & debugMode == F) {
          incProgress(0.2, detail = "Preparing to process original database")
          tmp <- tempfile()
          download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz.md5", tmp)
          clvMd5Serv <- readLines(tmp)
          saveRDS(clvMd5Serv, paste0("clvmd5_", format(Sys.time(), "%Y%m%d_%H%M%S"),".rds"))
          tmp <- tempfile()
          download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz", tmp)
          incProgress(0.4, detail = "Post-processing")
          clinvar <- as.data.table(read.table(gzfile(tmp)))
          colnames(clinvar) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
          clinvar$QUAL <- NULL
          clinvar$FILTER <- NULL
          clinvar$CHROM <- paste0("chr",clinvar$CHROM)
          clinvar$clnsig <- str_match(clinvar$INFO, "CLNSIG=\\s*(.*?)\\s*;")[,2]
          clinvar$clnrevstat <- str_match(clinvar$INFO, "CLNREVSTAT=\\s*(.*?)\\s*;")[,2]
          clinvar <- clinvar[!is.na(clinvar$clnsig),]
          clinvar$gene <- str_match(clinvar$INFO, "GENEINFO=\\s*(.*?)\\s*;")[,2]
          clinvar$INFO <- NULL
          clinvar <- clinvar %>% mutate(gene = str_split(clinvar$gene, "\\|")) %>%
            unnest(gene)
          clinvar$gene <- remove_gene_number(clinvar$gene)
          clinvar <- clinvar[,c("gene","CHROM","POS","ID","REF","ALT","clnsig","clnrevstat")]
          incProgress(0.1, detail = "Saving")
          sqldb_clinvar <- dbConnect(SQLite(), dbname=paste0("sqldb_clinvar_", format(Sys.time(), "%Y%m%d_%H%M%S")))
          dbWriteTable(sqldb_clinvar, "clinvar", clinvar, row.names=F, overwrite=T, append=F, field.types=NULL)
          clv_path <<- list.files(pattern = "sqldb_clinvar_", full.names = T)[length(list.files(pattern = "sqldb_clinvar_"))]
        }
      })
    }
    
    loadClinVar <- function(updateDb) {
        if (exists("gr_clinvar") == F) {
          withProgress(message = "Loading ClinVar database", {
            prepClinvar(updateDb = F)
          incProgress(0.5, detail = "Loading")
          sqldb_clinvar <- dbConnect(SQLite(), dbname=clv_path)
          gr_clinvar <<- GRanges(dbGetQuery(sqldb_clinvar, paste0('select CHROM from clinvar'))[,1],
                                 IRanges(dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar'))[,1],
                                         dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar'))[,1]+
                                           (nchar(dbGetQuery(sqldb_clinvar, paste0('select REF from clinvar'))[,1])-1)))
        })
      }
    }
    
    # COSMIC loader
    prepCosmic <- function() {
        # first, check if db exists or forced update, if not, create
        if (length(cmc_path) == 0) {
          withProgress(message = "Preparing COSMIC database", {
          # need RefSeq ex_by_ge for this
          loadRefSeq(force = F)
          incProgress(0.1, detail = "Preparing to process original database")
          library(Seurat)
          cmc_ori <- fread("db_ori/cmc_export.tsv")
          cmc_genes <- data.frame("ori" = setdiff(unique(cmc_ori$GENE_NAME), names(ex_by_ge)))
          cmc_genes$new <- cmc_genes$ori
          for (i in 1:nrow(cmc_genes)) {
            incProgress((0.6/nrow(cmc_genes)), detail = paste0("Update names: ", cmc_genes$ori[i]))
            cmc_genes$new[i] <- UpdateSymbolList(cmc_genes$ori[i])
          }
          incProgress(0.1, detail = "Post-processing")
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
            "GENE_NAME", "GENOMIC_MUTATION_ID", "Mutation CDS","Mutation AA","chr","start","end","COSMIC_SAMPLE_POSRATE","Mutation Description CDS","Mutation Description AA","MUTATION_SIGNIFICANCE_TIER", "CGC_TIER"
          )]
          # save to sql
          incProgress(0.1, detail = "Saving")
          sqldb_cosmic <- dbConnect(SQLite(), dbname=paste0("sqldb_cosmic_", format(Sys.time(), "%Y%m%d_%H%M%S")))
          dbWriteTable(sqldb_cosmic, "cosmic", cmc_ori, row.names=F, overwrite=T, append=F, field.types=NULL)
          # read coords, then disconnect, save path for later
          cmc_path <<- list.files(pattern = "sqldb_cosmic_", full.names = T)[length(list.files(pattern = "sqldb_cosmic_"))]
        })
      }
    }
    
    loadCosmic <- function() {
      withProgress(message = "Loading COSMIC database", {
        prepCosmic()
        # check if ranges are loaded, if not, load
        incProgress(0.1, detail = "Loading")
        if (exists("gr_cmc") == F) {
          incProgress(0.5, detail = "loading cosmic database")
          sqldb_cosmic <- dbConnect(SQLite(), dbname=cmc_path)
          gr_cmc <<- GRanges(dbGetQuery(sqldb_cosmic, paste0('select chr from cosmic'))[,1],
                             IRanges(dbGetQuery(sqldb_cosmic, paste0('select start from cosmic'))[,1],
                                     dbGetQuery(sqldb_cosmic, paste0('select end from cosmic'))[,1]))
        }
      })
    }
        
    # load and list panel files - need to transition to proper database at some point
    panelfiles <- list.files(path = "panels", pattern = ".panel", full.names = T)
    # maybe try this later  
    # panelfiles_table <- data.frame(date = str_extract(panel panelfiles, "\\d{8}"), time = str_extract(panelfiles, "(\\d+)(?!.*\\d)"))
    panelfiles_table <- as.data.frame(str_split_fixed(panelfiles, "_",2)) %>% group_by(V1) %>% arrange(desc(V2)) %>% dplyr::slice(1:1)
    panelfiles_table$V3 <- paste(panelfiles_table$V1, panelfiles_table$V2, sep = "_")
    panelfiles <- panelfiles[panelfiles %in% panelfiles_table$V3]
    dbx <- lapply(panelfiles, readRDS)
    
    # get all gene names- write this properly at a later point, vector needs to exist for the app not to crash on reactive drop down generation
    all_genes <- vector()
    for (i in 1:length(dbx)) {
      all_genes <- c(all_genes, dbx[[i]][["panelTable"]]$gene)
    }
    all_genes <- data.frame("gene" = sort(unique(all_genes)))
    
    # Check if PanelCat version is identical across panels, then finalise dbx by naming items, else skip loading and warn
    if (length(unique(unlist(nameList(dbx, "panelCatVersion")))) == 1) {
      names(dbx) <- nameList(dbx, "panelName")
      # how to check if panels with same name are completely different, e.g. compare rows of original target file, TBD
      status_dbload_success <- T
    } else {
      status_dbload_success <- F
    }
    # check if the reference databases versions of all loaded panels are the same
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
    exdb_path <- list.files(pattern = "exdb_", full.names = T)[length(list.files(pattern = "exdb_"))]
    cmc_path <- list.files(pattern = "sqldb_cosmic_", full.names = T)[length(list.files(pattern = "sqldb_cosmic_"))]
    clv_path <- list.files(pattern = "sqldb_clinvar_", full.names = T)[length(list.files(pattern = "sqldb_clinvar_"))]
    clv_md5_path <- list.files(pattern = "clvmd5", full.names = T)[length(list.files(pattern = "clvmd5"))]
    
    # create translate dataframe for silly panelTable variables (find better solution)
    inputxyc <- data.frame(code = names(dbx[[1]][["panelTable"]]), display = c("gene",
                                                                               "coding bases, not covered",
                                                                               "coding bases, covered",
                                                                               "coding bases, masked",
                                                                               "coding bases, covered not masked",
                                                                               "coding bases, total known",
                                                                               "coding bases, covered not masked  %",
                                                                               "coding bases, covered %",
                                                                               "cosmic muts, not covered",
                                                                               "cosmic muts, covered",
                                                                               "cosmic muts, masked",
                                                                               "cosmic muts, covered not masked",
                                                                               "cosmic muts, total known",
                                                                               "cosmic muts, covered not masked  %",
                                                                               "cosmic muts, covered %",
                                                                               "clinvar muts, not covered",
                                                                               "clinvar muts, covered",
                                                                               "clinvar muts, masked",
                                                                               "clinvar muts, covered not masked",
                                                                               "clinvar muts, total known",
                                                                               "clinvar muts, covered not masked  %",
                                                                               "clinvar muts, covered %"
    ))
    
    # JS URL render string
    render <- c(
      "function(data, type, row){",
      "  if(type === 'display'){",
      "    var a = '<a href=\"' + row[14] + '\" target=\"_blank\">' + data + '</a>';",
      "    return a;",
      "  } else {",
      "    return data;",
      "  }",
      "}"
    )
  }
  # Define UI ---------------------------------------------------------------
  
  css <- "div.dataTables_wrapper  div.dataTables_filter {width: 100%; float: none; text-align: center;}"
  
  ui <- ui <- fluidPage(
    useShinyjs(),
    tags$head(tags$style(HTML(css,"pre { white-space: pre-wrap; word-break: keep-all; }"))),
    headerPanel('PanelCat - GrCh37 (Research use only)'),
    
    sidebarPanel(
      conditionalPanel(
        condition = "input.test=='Gene metrics, X/Y'",
        uiOutput("xpanelsel"),
        selectInput("xcol", "Variable X",inputxyc$display[2:nrow(inputxyc)],selected = "coding bases, covered not masked"),
        uiOutput("ypanelsel"),
        selectInput("ycol", "Variable Y",inputxyc$display[2:nrow(inputxyc)],selected = "coding bases, covered not masked"),
        checkboxInput("showTipXY", "Show tool tips", value = T),
        span(tags$b("Coding bases:"),"RefSeq coding bases per gene, across all transcripts",tags$br(),tags$b("ClinVar muts:"),"counts of likely pathogenic / pathogenic ClinVar variants",tags$br(),tags$b("COSMIC muts:"),"counts of COSMIC Tier 1-3 census mutations",tags$br(),tags$br(),tags$b("Covered:"),"within targeted regions, but prior to applying mask regions (if provided).",tags$br(),tags$b("Not covered:"),"not within panel target regions.",tags$br(),tags$b("Covered, not masked:"), "within targeted regions, after applying the mask regions (if provided).",tags$br(),tags$b("Masked:"),"within mask regions (if provided)",tags$br(),tags$br(),tags$b("%:"),"referring to the total of coding bases / variants / mutations in a specific gene")
      ),
      conditionalPanel(
        condition = "input.test=='Gene metrics, column'",
        uiOutput("barpanelsel"),
        selectInput('dtset','Dataset', c("RefSeq coding bases","ClinVar mutations","COSMIC mutations")),
        textInput("diff_filter","Filter by fold change between 1st and 2nd selected panels", value = 1),
        selectInput("vsort", "Sort by", c("gene","fold change","total","ratio","max_cov","max_covp")),
        radioButtons("dtsetrel","display",c("absolute","relative"), selected = "absolute"),
        checkboxInput("showTipCol", "Show tool tips", value = T),
        span(tags$b("Coding bases:"),"RefSeq coding bases per gene, total of all transcripts",tags$br(),tags$b("ClinVar muts:"),"counts of likely pathogenic / pathogenic ClinVar variants",tags$br(),tags$b("COSMIC muts:"),"counts of COSMIC Tier 1-3 census mutations",tags$br(), tags$b("Filter by fold change:"),"If more than one panel is selected, the FC for each gene is calculated as (maximum coverage / minimum coverage) across all panels. '1' means no difference, and only values >1 will be effective.",tags$br(),tags$b("Sort by:"),"'fold change' refers to the difference between panels described above. 'Ratio' will sort the genes by the ratio of coverage between the first and the second selected panel. 'max_covp' will sort the genes by the maximum relative coverage of coding bases or mutations in any of the selected panels.")
      ),
      conditionalPanel(
        condition = "input.test=='Gene metrics, search'",
        uiOutput("genepanelsel"),
        uiOutput("genesel"),
        selectInput('dtsetgenes','Dataset', c("RefSeq coding bases","ClinVar mutations","COSMIC mutations")),
        radioButtons("dtsetrelg","display",c("absolute","relative"), selected = "absolute"),
        checkboxInput("showTipSearch", "Show tool tips", value = T),
        span(tags$b("Coding bases:"),"RefSeq coding bases per gene, total of all transcripts",tags$br(),tags$b("ClinVar muts:"),"counts of likely pathogenic / pathogenic ClinVar variants",tags$br(),tags$b("COSMIC muts:"),"counts of COSMIC Tier 1-3 census mutations")
      ),
      conditionalPanel(
        condition = "input.test == 'Exon table'",
        uiOutput("expanelsel"),
        downloadButton("exon_data_saved", "Download selected data"),
        checkboxInput("showTipExTab", "Show tool tips", value = T)
      ),
      conditionalPanel(
        condition = "input.test == 'Exon graph'",
        uiOutput("exoncomp_panelsel"),
        uiOutput("exoncovgsel"),
        uiOutput("exoncovtsel"),
        radioButtons("exoncomp_rel","display",c("absolute","relative"), selected = "absolute"),
        checkboxInput("showTipExGr", "Show tool tips", value = T)
      ),
      conditionalPanel(
        condition = "input.test == 'COSMIC table'",
        uiOutput("mutpanelsel"),
        checkboxInput("hideCmcBl","Hide masked", value = T),
        checkboxInput("showTipCosTab", "Show tool tips", value = T),
        span(tags$b("CDS: "),"Coding sequence",tags$br(), tags$br(), tags$b("AA: "),"Amino acid",tags$br(), tags$br(), tags$b("CMC Tier: "),"COSMIC cancer mutation census Tier",tags$br(), tags$br(), tags$b("CGC Tier: "),"COSMIC cancer gene census Tier")
      ),
      conditionalPanel(
        condition = "input.test == 'ClinVar table'",
        uiOutput("clvpanelsel"),
        checkboxInput("hideClvBl","Hide masked", value = T),
        checkboxInput("showTipClvTab", "Show tool tips", value = T),
        span(tags$b("REF: "),"Reference allele",tags$br(), tags$br(), tags$b("ALT: "), "Alternate allele", tags$br(), tags$br(), tags$b("ClnSig: "), "Clinical significance", tags$br(), tags$br(), tags$b("ClnRevStat: "),"Clinical review status")
      ),
      conditionalPanel(
        condition = "input.test=='Non-covered mutation rate'",
        selectInput("ncovSort","Sort by",c("panel","Target genes, with mask","Target genes, without mask","All genes, with mask","All genes, without mask"), selected = "panel"),
        checkboxInput("showTipNonCov", "Show tool tips", value = T),
        span(tags$b("Target genes, without mask:")," Considering only genes that are actually targeted by a panel, and without applying masking regions",tags$br(),tags$br(),tags$b("Target genes, with mask:"),"Considering only genes that are actually targeted by a panel, after applying masking regions (if provided)",tags$br(),tags$br(),tags$b("All genes, without mask: "),"Considering all genes in the COSMIC CMC database (including those that are not target genes of a panel), and without applying masking regions.",tags$br(),tags$br(),tags$b("All genes, with masking: "),"Considering all genes in the COSMIC CMC database (including those that are not target genes of a panel), after applying masking regions (if provided.")
      ),
      conditionalPanel(
        condition = "input.test=='Gene metrics, table'",
        uiOutput("sel_tablePanels"),
        selectInput('tableVars', 'Variables (choose multiple)', inputxyc$display[2:nrow(inputxyc)], selected = c("coding bases, covered not masked"), multiple = T),
        verbatimTextOutput("tableInfo2"),
        downloadButton("data_saved","Download selected data"),
        checkboxInput("showTipTable", "Show tool tips", value = T),
        span(tags$b("Coding bases:"),"RefSeq coding bases per gene, across all transcripts",tags$br(),tags$b("ClinVar muts:"),"counts of likely pathogenic / pathogenic ClinVar variants",tags$br(),tags$b("COSMIC muts:"),"counts of COSMIC Tier 1-3 census mutations",tags$br(),tags$br(),tags$b("Covered:"),"within targeted regions, but prior to applying mask regions (if provided).",tags$br(),tags$b("Not covered:"),"not within panel target regions.",tags$br(),tags$b("Covered, not masked:"), "within targeted regions, after applying the mask regions (if provided).",tags$br(),tags$b("Masked:"),"within mask regions (if provided)",tags$br(),tags$br(),tags$b("%:"),"referring to the total of coding bases / variants / mutations in a specific gene")
      ),
      conditionalPanel(
        condition = "input.test=='Panels'",
        uiOutput("sel_panelInspect")),
      conditionalPanel(
        condition = "input.test=='New analysis'",
        fileInput("file1","1: Panel target file (tab delimited, bed file)"),
        checkboxInput("zeroIndex","2: Zero-Indexed (.bed convention)?", value = T),
        textInput("pName","3: Panel Name (abbreviation)"),
        textInput("pfName","4: Panel Name (full)"),
        textInput("pRows","5: start row, end row"),
        textInput("pCols","6: Column order: Chr, Start, Stop"),
        fileInput("blackl","7: Panel mask file"),
        textInput("pmRows","8: Mask start row, end row"),
        textInput("pmCols","9: Mask column order: Chr, Start, Stop"),
        actionButton("startb", "10: Start!"),
        downloadButton("save_state", "11: Save to file"),
        fileInput("panelUp","12: processed panel file", accept = ".panel"),
        actionButton("panelUpButton","12: upload"),
        checkboxInput("showTipAnalysis", "Show tool tips", value = T)
        # comment the next line if hosting for others
        ,actionButton("updateb", "13: UPDATE ALL")
      ),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        id = "test",
        tabPanel("Gene metrics, X/Y",
                 htmlOutput("textXY"),
                 plotlyOutput("plot1"),
                 verbatimTextOutput("stats")),
        tabPanel("Gene metrics, column", 
                 htmlOutput("textCol"),
                 plotOutput("plot2")),
        tabPanel("Gene metrics, search", 
                 htmlOutput("textSearch"), 
                 verbatimTextOutput("unpanels"),
                 plotlyOutput("plotGenes")),
        tabPanel("Gene metrics, table", 
                 htmlOutput("textTable"),
                 DT::dataTableOutput("table", width = 100)),
        tabPanel("Exon graph",
                 htmlOutput("textExGra"),
                 plotlyOutput("exonCompP")),
        tabPanel("Exon table",
                 htmlOutput("textExTab"),
                 DT::dataTableOutput("exons", width = 100)),
        tabPanel("ClinVar table",
                 htmlOutput("textClvTab"),
                 DT::dataTableOutput("clvmuts", width = 100)),
        tabPanel("COSMIC table", 
                 htmlOutput("textCosTab"),
                 DT::dataTableOutput("cmcmuts", width = 100)),
        tabPanel("Non-covered mutation rate",
                 htmlOutput("textNonCov"),
                 plotlyOutput("nCovPosRate")),
        tabPanel("Panels", span("Inspect previously analyzed panels, including time stamps (i.e. versions) of the reference databases."),
                 verbatimTextOutput("viewP")),
        tabPanel("New analysis",
                 htmlOutput("textAnalysis"),
                 tableOutput("contents"),
                 tableOutput("maskFileContent")),
        tabPanel("Info",
                 span("PanelCat is an open-soure application designed to analyse and visualise NGS panel target regions, and store the analyses for quick access.", tagList(a("See GitHub repository.", href = "https://github.com/aoszwald/panelcat", target="blank")), tags$br(), "Research use only. PanelCat is released under AGPL-3 License. Author: André Oszwald."),
                 imageOutput("catpic"))
      ))
  )
  
  
  # # Define MESSY MESSYserver logic  --------------------------------------------------
  
  server <- function(input, output) {
    options(shiny.maxRequestSize=500*1024^2)
    
    observe({
      shinyalert("DISCLAIMER and user agreement:", text = "This software is intended for reasearch use only, and not intended to make medical decisions. By proceeding, the user agrees to take sole responsibility, and not to hold the authors/providers of this software responsible, for any decisions based on information obtained through this application.", confirmButtonText = "I understand and agree.", confirmButtonCol = "darkblue", closeOnEsc = F)
    })
    
    # check if DBs loaded
    if (status_dbload_success != T) {
      shinyalert(title = "Error loading panels", test = "Reason: Not all panels analyzed with same Version of PanelCat. Please restore latest functional state, update all, and then add new panels.")
    }
    
    if (status_dbload_conflict == T) {
      shinyalert(title = "Warning!", text = "Not all loaded panels were analyzed with identical database versions. Suggest to remove and re-analyze affected panels, or update all.")
    }
    
    if (length(cmc_path) == 0 & length(list.files(path = "db_ori", pattern="^cmc_export\\.tsv$")) == 0) {
      shinyalert(title = "Warning!", text  = "Please download COSMIC CMC database ('cmc_export.tsv' from https://cancer.sanger.ac.uk/cosmic) into the panelcat subdirectory 'db_ori'. You will need to create a COSMIC account. Proceeding without this file will lead to unexpected app behaviour.")
    }
    
    # define reactive input choices
    dbxn <- reactiveValues(panelNames = names(dbx), geneNames = all_genes)
    
    # tooltips
    output$textXY <- renderText({if (input$showTipXY == T) {HTML("Compare target region metrics of panels in an X/Y point graph (scatterplot). Each datapoint represents a gene and the corresponding coverage of total coding (exonic) bases, or known mutations from ClinVar or COSMIC databases. You can compare metrics between two different panels, or different metrics within the same panel.</br>Below the graph, you will find a text box contrasting the different sets of target genes.</br>The RefSeq metrics refer to the sum of all exon-coding bases of all transcripts of a gene. The ClinVar or COSMIC metrics include only pathogenic/likely pathogenic ClinVar variants and Tier 1-3 COSMIC census mutations (unlike the tabs 'ClinVar table' and 'COSMIC table').")} else {""}})
    output$textCol <- renderText({if (input$showTipCol == T) {HTML("Visualize coverage of coding bases, or Clinvar/COSMIC mutations all target genes of one or multiple panels and the extent of variant masking. Each horizontal column (or set of columns, if multiple panels are selected) represents one target gene.</br>The light gray bar represents the sum of all known exon bases of all transcripts, or the total number of known ClinVar / COSMIC mutations, per gene (i.e., the maximum that can be attained). The opaque bar in the foreground represents the coverage of coding bases or mutations by a selected panel. <b>The color mapping will change if panel selection is changed.</b></br>Masked bases or variants are indicated by a lighter shade of same color of the foreground. If you are unsure what this looks like, select the TSO500 and look at the genes KMT2B-D.</br>The ClinVar or COSMIC metrics include only pathogenic/likely pathogenic ClinVar variants and Tier 1-3 COSMIC census mutations (unlike the tabs 'ClinVar table' and 'COSMIC table').</br>You can right-click the graph and save image as a .png image file. Do not select too many panels, or performance will suffer. Consider using the 'Gene Metrics, Search' function instead.")} else {""}})
    output$textSearch <- renderText({if (input$showTipSearch == T) {HTML("Visualize coverage of coding bases, or Clinvar/COSMIC mutations of one or more genes of interest across one or multiple panels. Each horizontal column (or set of columns, if multiple panels are selected) represents one target gene. </br>The light gray shaded bar represents the sum of all known exon bases of all transcripts, or the total number of known ClinVar / COSMIC mutations, per gene (i.e., the maximum that can be attained). The opaque bar in the foreground represents the coverage by a selected panel. <b>The color mapping will change if panel selection is changed.</b></br>Masked bases or variants are indicated by a lighter shade of same color of the foreground. If you are unsure what this looks like, select the TSO500 and look at the genes KMT2B-D</br>The ClinVar or COSMIC metrics include only pathogenic/likely pathogenic ClinVar variants and Tier 1-3 COSMIC census mutations (unlike the tabs 'ClinVar table' and 'COSMIC table').</br>Since only selected genes are displayed at a time, performance is acceptable when comparing large numbers of panels.")} else {""}})
    output$textTable <- renderText({if (input$showTipTable == T) {HTML("View the data used to construct the graphs in tabular form. You can select as many panels and as many metrics as you like to compare (but beware, the table will eventually become very wide).</br>The ClinVar or COSMIC metrics include only pathogenic/likely pathogenic ClinVar variants and Tier 1-3 COSMIC census mutations (unlike the tabs 'ClinVar table' and 'COSMIC table').</br>You may download the currently selected data in tabular form (the export will be filtered based on your input).")} else {""}})
    output$textNonCov <- renderText({if (input$showTipNonCov == T) {HTML("View the estimated rate of tested samples harboring COSMIC Tier 1-3 census mutations that will not be detected with a panel because they lie outside the specified target regions. The estimation is derived from the number of positive samples, and the number of samples tested for this mutation, documented in the COSMIC database. The estimate does not account for different tumor entities, or whether samples were analysed using either genome-wide or targeted screens.")} else {""}})
    output$textExGra <- renderText({if (input$showTipExGr == T) {HTML("Visualise coverage of individual exons of any transcript in one or multiple panels.</br>The light gray shaded bar represents the sum of all known exon bases of all transcripts, or the total number of known ClinVar / COSMIC mutations, per gene (i.e., the maximum that can be attained). The opaque bar in the foreground represents the coverage by a selected panel.</br>Masked bases or variants are indicated by a lighter shade of same color of the foreground. If you are unsure what this looks like, select the TSO500 and look at the genes KMT2B-D.</br> IMPORTANT HINT: select the panels you wish to compare FIRST, because selecting new panels will discard your current gene/transcript selection.")} else {""}})
    output$textExTab <- renderText({if (input$showTipExTab == T) {HTML("For a specific panel, assess the targeted portion of all exons of all transcripts of all genes listed in refseq. Exon coverage % is rounded to five decimals.</br>Hint: Try filtering the columns for specific genes and transcripts, and searching the table for specific mutations. Clicking the transcript accession number will open the corresponding NCBI entry in your browser. You may download the currently selected data in tabular form (the export will be filtered based on your input).")} else {""}})
    output$textCosTab <- renderText({if (input$showTipCosTab == T) {HTML("Inspect coding mutations (from the COSMIC mutation census database) that are targeted by a specific panel. This table includes ALL coding mutations (i.e. not only Tier 1-3 census mutations), unlike the other tabs in PanelCat. The COSMIC frequency is rounded to two significant digits for each entry.</br>Hint: Try filtering the columns for specific genes and transcripts, and searching the table for specific mutations. Clicking the ID will open the corresponding COSMIC entry in your browser.</br>The COSMIC frequency is derived from the number of positive samples, and the number of samples tested for this mutation, documented in the COSMIC database. The estimate does not account for different tumor entities, or whether samples were analysed using either genome-wide or targeted screens.")} else {""}})
    output$textClvTab <- renderText({if (input$showTipClvTab == T) {HTML("Inspect ClinVar mutations that are targeted by a specific panel. This table includes ALL mutations (i.e. not only pathogenic / likely pathogenic), unlike the other tabs in PanelCat.</br>Hint: Try filtering the columns for specific genes and transcripts, and searching the table for specific mutations. Clicking the ID will open the corresponding ClinVar entry in your browser.")} else {""}})
    output$textAnalysis <- renderText({if (input$showTipAnalysis == T) {HTML("To create a new panel file, choose a tab-delimited target region file (1). These files are typically provided as .bed files by the manufacturers of NGS panels.</br> 
                                                                             If the target region file does not adhere to .bed conventions <a href = 'https://en.wikipedia.org/wiki/BED_(file_format)#Coordinate_system' target='blank'>(read this)</a>, uncheck the box (2). One indication of this case is if your target region file is a .txt file (and not a .bed file)</br>
                                                                             Specify an abbreviation for use in drop-down menus, charts, etc. (3).</br>
                                                                             Specify a full name for a meaningful identification (4).</br>
                                                                             If the file contains additional headers or other rows at the end, specify the line numbers (e.g. 2,2000) of file to be included in the analysis (5).</br>
                                                                             If the first three columns of the file do not correspond to chromosome, start and stop coordinates, specifiy the position and order of the respective columns (e.g. 2,3,4) to be included in the analysis (6).</br>
                                                                             Use the same approach for an optional mask file (7-9).</br>
                                                                             For convenience, panel analyses can be stored by download (11) after analysis is complete and uploaded again at a later time (12).</br>
                                                                             You can initiate an update of RefSeq and ClinVar databases and a re-analysis of all currently loaded panels (13).</br>
                                                                             To update COSMIC, you will have to remove the existing sqldb_cosmic_(date)_(time) file and place a current cmc_export.tsv in the db_ori folder.")} else {""}})
  
      # define reactive inputs
    output$xpanelsel <- renderUI({
      selectInput('xpanel', 'Panel X', dbxn$panelNames, selected = dbxn$panelNames[1])
    })
    
    output$ypanelsel <- renderUI({
      selectInput('ypanel', 'Panel Y', dbxn$panelNames, selected = dbxn$panelNames[2])
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
    
    output$exoncomp_panelsel <- renderUI({
      selectInput('exoncomp_panel', 'Panels (choose multiple)', dbxn$panelNames, selected = dbxn$panelNames[1], multiple = T)
    })
    output$exoncovgsel <- renderUI({
      validate(need(input$exoncomp_panel, "Select at least one panel"))
      selectInput('exoncovg', 'Gene', exoncovd()[["group_name"]], selected = exoncovd()[["group_name"]][1])
    })
    output$exoncovtsel <- renderUI({
      validate(need(input$exoncovg,"Select a gene"))
      selectInput('exoncovt', 'Transcript', exoncovd1()[["transcript"]], selected = exoncovd1()[["transcript"]][1])
    })
    output$sel_panelInspect <- renderUI({
      selectInput('panelInspect', "Select Panel", dbxn$panelNames, selected = dbxn$panelNames[1], size = 25, selectize = F)
    })
    
    output$sel_tablePanels <- renderUI({
      selectInput('tablePanels', 'Panels (choose multiple)', dbxn$panelNames, selected = dbxn$panelNames[1], multiple = T)
    })
    
    output$catpic <- renderImage(list(src = "cats.png", width = 600, height = 400), deleteFile = F)
    
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
        if (panelUpFile()[["panelName"]] %in% names(dbx)) {
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
          
          dbxn$panelNames <- c(dbxn$panelNames, panelUpFile()[["panelName"]])
          dbxn$geneNames <- all_genes
          shinyalert(title = "complete")
        }
      }
    })
    
    
    
    # X/Y ---------------------------------------------------------------------
    
    # scatter data
    scatter_data <- reactive({
      validate(need(input$xpanel, 'Please select panel X'),need(input$ypanel, 'Please select panel Y'))
      temp <- rbindlist(sapply(dbx[c(input$xpanel, input$ypanel)], "[", "panelTable"), idcol = "panel") %>%
        mutate(panel = gsub(".panelTable","",panel))
      
      temp <- full_join(data.frame(gene = temp[panel == input$xpanel][["gene"]],
                                   x = temp[panel == input$xpanel][[inputxyc$code[match(input$xcol, inputxyc$display)]]]),
                        data.frame(gene = temp[panel == input$ypanel][["gene"]],
                                   y = temp[panel == input$ypanel][[inputxyc$code[match(input$ycol, inputxyc$display)]]])) %>% replace(is.na(.),0)
      temp
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
    
    output$stats <- renderText({
      validate(need(input$xpanel, "Please select panel X"), need(input$ypanel, "Please select panel Y"))
      paste0(input$xpanel,": ",nrow(as.data.frame(dbx[[input$xpanel]][["panelTable"]]))," genes \n",
                                      input$ypanel,": ",nrow(as.data.frame(dbx[[input$ypanel]][["panelTable"]]))," genes \n\n",
                                      "Genes exclusively in Panel ",input$xpanel,": ",setAB(),"\n\n","Genes exclusively in Panel ",input$ypanel,": ",setBA(),"\n\n","Genes shared: ",isectBA())})
    
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
      validate(
        need(input$panelbar, 'Please select at least one panel!'),
      )
      bars_cov1 <- rbindlist(sapply(dbx[input$panelbar], "[", "panelTable"), idcol = "panel") %>%
        filter (case_when(input$dtset == "ClinVar mutations" ~ clv_tot != 0,
                          input$dtset == "RefSeq coding bases" ~ !is.na(pcb_tot),
                          input$dtset == "COSMIC mutations" ~ cmc_tot != 0)) %>% 
        as_tibble() %>% 
        mutate_at(c("panel", "gene"), as.factor) %>% 
        complete(panel,gene) %>%
        mutate(panel = gsub(".panelTable","",panel))
      bars_cov1 <- left_join(bars_cov1, bars_cov1 %>%
                               replace(is.na(.), 0) %>%
                               group_by(gene)%>%
                               summarise(FC = case_when(input$dtset == "RefSeq coding bases" ~ max(pcb_covp) / min(pcb_covp),
                                                           input$dtset == "ClinVar mutations" ~ max(clv_covp) / min(clv_covp),
                                                           input$dtset == "COSMIC mutations" ~ max(cmc_covp) / min(cmc_covp)),
                                         total = case_when(input$dtset == "RefSeq coding bases" ~ max(pcb_tot),
                                                           input$dtset == "ClinVar mutations" ~ max(clv_tot),
                                                           input$dtset == "COSMIC mutations" ~ max(cmc_tot)),
                                         ratio = case_when(input$dtset == "RefSeq coding bases" ~ pcb_covp[panel == input$panelbar[1]]/pcb_covp[panel == input$panelbar[2]],
                                                           input$dtset == "ClinVar mutations" ~ clv_covp[panel == input$panelbar[1]]/clv_covp[panel == input$panelbar[2]],
                                                           input$dtset == "COSMIC mutations" ~ cmc_covp[panel == input$panelbar[1]]/cmc_covp[panel == input$panelbar[2]]),
                                         max_cov = case_when(input$dtset == "RefSeq coding bases" ~ max(pcb_cov),
                                                             input$dtset == "ClinVar mutations" ~ max(clv_cov),
                                                             input$dtset == "COSMIC mutations" ~ max(cmc_cov)),
                                         max_covp = case_when(input$dtset == "RefSeq coding bases" ~ max(pcb_covp),
                                                              input$dtset == "ClinVar mutations" ~ max(clv_covp),
                                                              input$dtset == "COSMIC mutations" ~ max(cmc_covp)))) %>%
        filter(FC >= as.numeric(input$diff_filter) | is.nan(FC))
      
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
        summarise(width = max(case_when(input$dtset == "RefSeq coding bases" ~ pcb_tot,
                                        input$dtset == "ClinVar mutations" ~ clv_tot,
                                        input$dtset == "COSMIC mutations" ~ cmc_tot)))
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
                                                             input$dtsetrel == "relative" ~100),
                                     fill = "total")) +
        geom_col(d = bars_cov(), aes(x = gene, y = case_when(input$dtset == "RefSeq coding bases" & input$dtsetrel == "absolute" ~ pcb_covt,
                                                             input$dtset == "ClinVar mutations" & input$dtsetrel == "absolute" ~ clv_covt,
                                                             input$dtset == "COSMIC mutations" & input$dtsetrel == "absolute" ~ cmc_covt,
                                                             input$dtset == "RefSeq coding bases" & input$dtsetrel == "relative" ~ pcb_covtp,
                                                             input$dtset == "ClinVar mutations" & input$dtsetrel == "relative" ~ clv_covtp,
                                                             input$dtset == "COSMIC mutations" & input$dtsetrel == "relative" ~ cmc_covtp)
                                     , fill = panel), position = "dodge", alpha = 0.5) +
        geom_col(d = bars_cov(), aes(x = gene,  y = case_when(input$dtset == "RefSeq coding bases" & input$dtsetrel == "absolute" ~ pcb_cov,
                                                              input$dtset == "ClinVar mutations" & input$dtsetrel == "absolute" ~ clv_cov,
                                                              input$dtset == "COSMIC mutations" & input$dtsetrel == "absolute" ~ cmc_cov,
                                                              input$dtset == "RefSeq coding bases" & input$dtsetrel == "relative" ~ pcb_covp,
                                                              input$dtset == "ClinVar mutations" & input$dtsetrel == "relative" ~ clv_covp,
                                                              input$dtset == "COSMIC mutations" & input$dtsetrel == "relative" ~ cmc_covp), 
                                     fill = panel), color = "black", size = 0.1, position = "dodge") +
        ylab(case_when(input$dtset == "RefSeq coding bases" & input$dtsetrel == "absolute" ~ "protein coding bases, absolute",
                       input$dtset == "ClinVar mutations" & input$dtsetrel == "absolute" ~ "ClinVar variants, absolute",
                       input$dtset == "COSMIC mutations" & input$dtsetrel == "absolute" ~ "COSMIC mutations, absolute",
                       input$dtset == "RefSeq coding bases" & input$dtsetrel == "relative" ~ "protein coding bases, relative",
                       input$dtset == "ClinVar mutations" & input$dtsetrel == "relative" ~ "ClinVar variants, relative",
                       input$dtset == "COSMIC mutations" & input$dtsetrel == "relative" ~ "COSMIC mutations, relative")) +
        theme_minimal(base_size = 15) +
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
      rbindlist(sapply(dbx[input$panelgenes], "[", "panelTable"), idcol = "panel") %>%
        filter(gene %in% input$genes) %>%
        mutate(panel = gsub(".panelTable","",panel)) %>%
        filter (case_when(input$dtset == "ClinVar mutations" ~ clv_tot != 0,
                          input$dtset == "RefSeq coding bases" ~ !is.na(pcb_tot),
                          input$dtset == "COSMIC mutations" ~ cmc_tot != 0)) %>% 
        as_tibble() %>% 
        mutate_at(c("panel", "gene"), as.factor) %>% 
        complete(panel,gene) %>%
        transform(panel = reorder(.$panel, rev(order(.$panel))))
    })
    
    genes_tot <- reactive({
      bars_gene() %>%
        replace(is.na(.), 0) %>%
        group_by(gene)%>%
        summarise(width = max(case_when(input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "absolute" ~ pcb_tot,
                                        input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "absolute" ~ clv_tot,
                                        input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "absolute" ~ cmc_tot,
                                        input$dtsetrelg == "relative" ~100)))
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
      genegraph <- ggplotly(
        ggplot() +
          geom_col(d = genes_tot(), aes(y = gene, x = width, fill = "total")) +
          geom_col(d = bars_gene(), aes(y = gene, x = case_when(input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "absolute" ~ pcb_covt,
                                                                input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "absolute" ~ clv_covt,
                                                                input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "absolute" ~ cmc_covt,
                                                                input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "relative" ~ pcb_covtp,
                                                                input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "relative" ~ clv_covtp,
                                                                input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "relative" ~ cmc_covtp), 
                                        fill = panel), position = "dodge", alpha = 0.3) +
          geom_col(d = bars_gene(), aes(y = gene,  x = case_when(input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "absolute" ~ pcb_cov,
                                                                 input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "absolute" ~ clv_cov,
                                                                 input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "absolute" ~ cmc_cov,
                                                                 input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "relative" ~ pcb_covp,
                                                                 input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "relative" ~ clv_covp,
                                                                 input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "relative" ~ cmc_covp), 
                                        fill = panel), color = "black", size = 0.05, position = "dodge") +
          xlab(case_when(input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "absolute" ~ "protein coding bases, absolute",
                         input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "absolute" ~ "ClinVar variants, absolute",
                         input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "absolute" ~ "COSMIC mutations, absolute",
                         input$dtsetgenes == "RefSeq coding bases" & input$dtsetrelg == "relative" ~ "protein coding bases, relative",
                         input$dtsetgenes == "ClinVar mutations" & input$dtsetrelg == "relative" ~ "ClinVar variants, relative",
                         input$dtsetgenes == "COSMIC mutations" & input$dtsetrelg == "relative" ~ "COSMIC mutations, relative")) +
          theme_minimal(base_size = 15) +
          ggtitle("Coverage (opaque) / masked (shaded)") +
          scale_fill_manual(values = genes_colvec()),
        height = case_when((10 * length(input$genes) * length(input$panelgenes)) <= 300 ~ 300,
                           (10 * length(input$genes) * length(input$panelgenes)) > 300 ~ (10 * length(input$genes) * length(input$panelgenes)))
      ) %>% remove_alpha_legend_genes_reverse()
    })
    
    unpanelsel <- reactive({
      bars_gene_full <- bars_gene()
      bars_gene_Na <- unique(bars_gene_full$panel[is.na(bars_gene_full$pcb_covt)])
      setdiff(unique(bars_gene_full$panel), bars_gene_Na)
    })
    
    # union panels text
    output$unpanels <-  renderText(c("Panels containing all selected genes: ", paste(unpanelsel(), collapse = " ")))
    
    
    # Non-targeted mutation stats -----------------------------------------------
    
    nCovPosRate_df <- reactive({
      nCovPosRate <- data.table("panel" = unlist(nameList(dbx, "panelName")),
                                "All genes, with mask" = unlist(nameList(dbx, "cmcNcovPosRateTotal_bl")),
                                "All genes, without mask" = unlist(nameList(dbx, "cmcNcovPosRateTotal")),
                                "Target genes, with mask" = unlist(nameList(dbx, "cmcNcovPosRate_bl")), 
                                "Target genes, without mask" = unlist(nameList(dbx, "cmcNcovPosRate")))
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
        ggplot(d = nCovPosRate_df(), aes(y = panel, x = value*100, fill = variable)) +
          geom_col(position = position_dodge2(preserve = "single")) +
          labs(x = paste0("Estimated number of Tier 1-3 COSMIC census mutations", "\n","outside of panel target regions, per 100 samples tested")) +
          theme_minimal(base_size = 12),
        height = (10*nrow(nCovPosRate_df()) + 200)
      ) %>% reverse_legend_labels()
    })
      
    # Exon coverage table ----------------------------------------------------- 
    
    table_exons <- reactive({
      loadExDb()
      left_join(dbx[[input$expanel]][["exon_coverage"]], ex_by_ge_df1) %>% 
        mutate(covp = cov_width / width, covtp = covt_width / width, group_name = as.factor(group_name), transcript = as.factor(transcript), href = paste("https://www.ncbi.nlm.nih.gov/nuccore/",transcript)) %>%
        setNames(c("gene","chromosome","start","end","strand","width","covered bases","covered bases, excluding masked","exon_id","transcript","exon","covered bases %, excluding masked","covered bases %", "href"))
        })
    
    output$exons <- DT::renderDataTable({
      DT::datatable(table_exons(), escape = F,
                    filter = list(position = "top", clear = F),
                    options = list(search = list(regex = TRUE, caseInsensitive = T), iDisplayLength = 100,
                                   columnDefs = list(
                                     list(targets = 10, render = DT::JS(render)),
                                     list(targets = 14, visible = F)
                                   ))) %>%
        formatRound(columns = c(12,13), digits = 5)
    })
    
    output$exon_data_saved = downloadHandler(
      filename = function(){
        paste0(input$expanel, "_exons_",format(Sys.time(), "%Y%m%d_%H%M%S"),".csv")},
      content = function(file){
        fwrite(table_exons()[input[["exons_rows_all"]], ], file)
      }
    )
    
    # Exon coverage graph -----------------------------------------------------
    
    exoncovd <- reactive({
      loadExDb()
      rbindlist(sapply(dbx[input$exoncomp_panel], "[", "exon_coverage"), idcol = "panel") %>%
      mutate(panel = gsub(".exon_coverage","",panel))
    })
    
    exoncovd1 <- reactive({
      exoncovd() %>% filter(group_name == input$exoncovg) %>%
        left_join(ex_by_ge_df1, relationship = "many-to-many")
    })
    
    exoncovd2 <- reactive({
      exoncovd1() %>% filter(transcript == input$exoncovt) %>%
        mutate(covp = (cov_width / width)*100, covtp = (covt_width / width)*100, exon = as.factor(exon)) %>%
        transform(exon = reorder(exon, rev(order(sort(exon)))))
    })
    
    exons_colvec <- reactive({
      setNames(c("grey95",distinctColorPalette(length(unique(exoncovd2()[["panel"]])))),
               c("total",unique(exoncovd2()[["panel"]])))
    })
    
    output$exonCompP <- renderPlotly({
      validate(
        need(input$exoncomp_panel, 'Please select at least one panel'),
        need(input$exoncovg, 'Please select one gene'),
        need(input$exoncovt, 'Please select one trx')
      )
      exonplot <- ggplotly(
        ggplot() +
          geom_col(d = exoncovd2() %>% 
                     group_by(exon) %>%
                     summarise(width = max(width)), aes(x = exon, y = case_when(input$exoncomp_rel == "absolute" ~ width,
                                                                                input$exoncomp_rel == "relative" ~100),
                                                        fill = "total")) +
          geom_col(d = exoncovd2(), aes(x = exon, y = case_when(input$exoncomp_rel == "absolute" ~ covt_width,
                                                                input$exoncomp_rel == "relative" ~ covtp),
                                        fill = panel), position = "dodge", alpha = 0.2, show.legend = F) +
          geom_col(d = exoncovd2(), aes(x = exon, y = case_when(input$exoncomp_rel == "absolute" ~ cov_width,
                                                                input$exoncomp_rel == "relative" ~ covp), 
                                        fill = panel), color = "black", size = 0.05, position = "dodge") +
          ylab(case_when(input$exoncomp_rel == "absolute" ~ "Exon coverage, basepairs",
                         input$exoncomp_rel == "relative" ~ "Exon coverage, % bp")) +
          theme_minimal() +
          coord_flip() +
          ggtitle("Coverage (opaque) / masked (shaded)") +
          scale_fill_manual(values = exons_colvec()),
        height = case_when((10 * length(unique(exoncovd2()[["exon"]])) * length(input$exoncomp_panel)) <= 300 ~ 400,
                           (10 * length(unique(exoncovd2()[["exon"]])) * length(input$exoncomp_panel)) > 300 ~ (12 * length(unique(exoncovd2()[["exon"]])) * length(input$exoncomp_panel)))
        
      ) %>% remove_alpha_legend(input$exoncomp_panel)
    })
    
    # COSMIC coverage table -------------------------------------------------
    
    table_muts <- reactive({
      loadCosmic()
      sqldb_cosmic <- dbConnect(SQLite(), dbname=cmc_path)
      gr_test <- GRanges(dbx[[input$mutpanel]][["panelBed_input"]][["V1"]],
                         IRanges(
                           dbx[[input$mutpanel]][["panelBed_input"]][["V2"]],
                           dbx[[input$mutpanel]][["panelBed_input"]][["V3"]]))
      if (input$hideCmcBl == T | length(dbx[[input$mutpanel]][["blacklist"]]) == 1) {
        muts_overlaps <- findOverlaps(gr_cmc, gr_test, type = "within")@from
        dbGetQuery(sqldb_cosmic, paste0('SELECT * FROM cosmic WHERE rowid IN (', paste(muts_overlaps, collapse = ","),')')) %>%
        mutate(GENOMIC_MUTATION_ID = sprintf('<a href="https://cancer.sanger.ac.uk/cosmic/search?q=%s" target="_blank"> %s </a>',GENOMIC_MUTATION_ID,GENOMIC_MUTATION_ID),
               `Mutation Description CDS` = as.factor(`Mutation Description CDS`), `Mutation Description AA` = as.factor(`Mutation Description AA`)) %>%
        setNames(c("gene","Mutation ID","Mutation CDS","Mutation AA","chrom","start","end","COSMIC frequency","Mutation Description CDS","Mutation Description AA","CMC Tier","CGC Tier"))
      } else {
        gr_bl <- GRanges(dbx[[input$mutpanel]][["blacklist"]][["V1"]],
                         IRanges(
                           dbx[[input$mutpanel]][["blacklist"]][["V2"]],
                           dbx[[input$mutpanel]][["blacklist"]][["V3"]]))
        gr_test_bl <- unlist(GenomicRanges::subtract(gr_test, gr_bl))
        muts_overlaps <- findOverlaps(gr_cmc, gr_test_bl, type = "within")@from
        dbGetQuery(sqldb_cosmic, paste0('SELECT * FROM cosmic WHERE rowid IN (', paste(muts_overlaps, collapse = ","),')')) %>%
        mutate(GENOMIC_MUTATION_ID = sprintf('<a href="https://cancer.sanger.ac.uk/cosmic/search?q=%s" target="_blank"> %s </a>',GENOMIC_MUTATION_ID,GENOMIC_MUTATION_ID),
        `Mutation Description CDS` = as.factor(`Mutation Description CDS`), `Mutation Description AA` = as.factor(`Mutation Description AA`)) %>%
        setNames(c("gene","Mutation ID","Mutation CDS","Mutation AA","chrom","start","end","COSMIC frequency","Mutation Description CDS","Mutation Description AA","CMC Tier","CGC Tier"))
      }
    })
    
    output$cmcmuts <- DT::renderDataTable({
      DT::datatable(table_muts(), escape = F,
                    filter = list(position = "top", clear = F),
                    rownames = F,
                    options = list(search = list(regex = TRUE, caseInsensitive = T), iDisplayLength = 100)) %>%
        formatSignif(columns = c('COSMIC frequency'), digits = 2)
    })
    
    # CLINVAR coverage table -------------------------------------------------
    
    
    table_clvmuts <- reactive({
      loadClinVar(updateDb = FALSE)
      sqldb_clinvar <- dbConnect(SQLite(), dbname=clv_path)
      gr_test <- GRanges(dbx[[input$clvpanel]][["panelBed_input"]][["V1"]],
                         IRanges(
                           dbx[[input$clvpanel]][["panelBed_input"]][["V2"]],
                           dbx[[input$clvpanel]][["panelBed_input"]][["V3"]]))
      if (input$hideClvBl == T | length(dbx[[input$clvpanel]][["blacklist"]]) == 1) {
        muts_overlaps <- findOverlaps(gr_clinvar, gr_test, type = "within")@from
        dbGetQuery(sqldb_clinvar, paste0('SELECT * FROM clinvar WHERE rowid IN (', paste(muts_overlaps, collapse = ","),')')) %>%
          mutate(gene = as.factor(gene), CHROM = as.factor(CHROM), clnsig = as.factor(clnsig), clnrevstat = as.factor(clnrevstat), ID = sprintf('<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term=%s" target="_blank"> %s </a>',ID,ID)) %>%
          setNames(c("gene","chrom","start","ID","REF","ALT","ClinSig","ClnRevStat"))
      } else {
        gr_bl <- GRanges(dbx[[input$clvpanel]][["blacklist"]][["V1"]],
                         IRanges(
                           dbx[[input$clvpanel]][["blacklist"]][["V2"]],
                           dbx[[input$clvpanel]][["blacklist"]][["V3"]]))
        gr_test_bl <- unlist(GenomicRanges::subtract(gr_test, gr_bl))
        muts_overlaps <- findOverlaps(gr_clinvar, gr_test_bl, type = "within")@from
        dbGetQuery(sqldb_clinvar, paste0('SELECT * FROM clinvar WHERE rowid IN (', paste(muts_overlaps, collapse = ","),')')) %>%
          mutate(gene = as.factor(gene), CHROM = as.factor(CHROM), clnsig = as.factor(clnsig), clnrevstat = as.factor(clnrevstat), ID = sprintf('<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term=%s" target="_blank"> %s </a>',ID,ID)) %>%
          setNames(c("gene","chrom","start","ID","REF","ALT","ClinSig","ClnRevStat"))
      }
    })
    
    output$clvmuts <- DT::renderDataTable({
      DT::datatable(table_clvmuts(), escape = F,
                    filter = list(position = "top", clear = F),
                    rownames = F,
                    options = list(search = list(regex = TRUE, caseInsensitive = T), iDisplayLength = 100))
    })
    
    
    # Panel coverage Tables ------------------------------------------------------------------
    
    table_output <- reactive({
      table <- dcast(rbindlist(sapply(dbx[input$tablePanels], "[", "panelTable"), idcol = "panel") %>%
                       mutate(panel = gsub(".panelTable","",panel)),
                     gene~panel, value.var = inputxyc$code[match(input$tableVars, inputxyc$display)])
      table <- table[apply(table, 1, count_nna_func) > 1,]
      tablenames <- expand.grid(input$tablePanels, input$tableVars)
      colnames(table) <- c("gene", paste(tablenames$Var1, tablenames$Var2))
      table
    })
    
    output$table <- DT::renderDataTable({
      validate(
        need(input$tablePanels, 'Select at least one panel'),
        need(input$tableVars, "Select at least one variable")
      )
      DT::datatable(table_output(), rownames = F, filter = list(position = "top", clear = F),
                    options = list(search = list(regex = TRUE, caseInsensitive = T), iDisplayLength = 100))
    })
    
    output$data_saved = downloadHandler(
      filename = function(){
        paste0(paste(input$tablePanels, collapse = "_"),"_metrics_",format(Sys.time(), "%Y%m%d_%H%M%S"),".csv")},
      content = function(file){
        fwrite(table_output()[input[["table_rows_all"]], ], file)
      }
    )
    
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
        colVecs <- seq(1,3)
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
      trget$V2 <- as.numeric(trget$V2)
      trget$V3 <- as.numeric(trget$V3)
      
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
      
      
      if (input$pName %in% names(dbx)) {
        shinyalert(title = "Panel Name already in use. Please provide a new unique name.")
      } else if (is.null(input$pName)) {
        shinyalert(title = "Please provide a new unique panel name.")
      } else if ((sum(anyNA(as.numeric(trget$V2)), anyNA(as.numeric(trget$V3))) != 0)) {
        shinyalert(title = "Specified panel file columns or rows contain non-numeric entries. Please review.")
      } else if (blacklist_check == F) {
        shinyalert(title = "Specified mask file columns or rows contain non-numeric entries. Please review.")
      } else {
        
        withProgress(message = "Processing", {
          # load databases
          loadRefSeq(force = F)
          prepClinvar(updateDb = F)
          prepCosmic()
          
          # make Granges of panel and mask
          gr_test <- reduce(GRanges(trget$V1, IRanges(trget$V2, trget$V3)))
          
          if (!is.null(panelInput$mask)) {
            blacklist$V2 <- as.numeric(blacklist$V2)
            blacklist$V3 <- as.numeric(blacklist$V3)
            names(blacklist) <- c("V1", "V2","V3")
            gr_blacklist <- GRanges(blacklist$V1, IRanges(blacklist$V2, blacklist$V3))
          } else {
            gr_blacklist <- GRanges()
          }
          
          gr_test_bl <- reduce(unlist(subtract(gr_test, gr_blacklist)))
          
          incProgress(0.1, detail = "Find targeted exons")
          
          ### RefSeq
          # find target genes
          gl_exp <- data.frame("refseq" = unique(names(unlist(ex_by_ge)[findOverlaps(gr_test, unlist(ex_by_ge))@to])))
          
          # identify unique exons and their coverage
          incProgress(0.1, detail = "Find unique exon coverage (optional)")
          exons <- ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq]
          ex_fo <- findOverlaps(gr_test, exons)
          ex_fo_bl <- findOverlaps(gr_test_bl, exons)
          ex_by_ge_all_df <- subset(
            left_join(as.data.table(exons), left_join(
              as.data.table(pintersect(gr_test[queryHits(ex_fo)], exons[subjectHits(ex_fo)])) %>% 
                filter(hit == T) %>%
                group_by(exon_id) %>%
                summarise(covt_width = sum(width)),
              as.data.table(pintersect(gr_test_bl[queryHits(ex_fo_bl)], exons[subjectHits(ex_fo_bl)])) %>% 
                filter(hit == T) %>%
                group_by(exon_id) %>%
                summarise(cov_width = sum(width))),
            by = "exon_id"),
            select = c(group_name, seqnames, start, end, strand, width, covt_width, cov_width)) %>%
            distinct() %>%
            replace(is.na(.), 0)
          
          # reduce exon ranges
          incProgress(0.1, detail = "Reducing exon GRanges")
          exons_red <- reduce(ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq])
          
          ex_fo1 <- findOverlaps(gr_test, exons_red)
          ex_fo_pint1 <- pintersect(gr_test[queryHits(ex_fo1)], exons_red[subjectHits(ex_fo1)])
          table_refseq_cov <- as.data.table(ex_fo_pint1) %>% group_by(group_name) %>% summarise(pcb_covt = sum(width))
          
          ex_fo2 <- findOverlaps(gr_test_bl, exons_red)
          ex_fo_pint2 <- pintersect(gr_test_bl[queryHits(ex_fo2)], exons_red[subjectHits(ex_fo2)])
          table_refseq_covbl <- as.data.table(ex_fo_pint2) %>% group_by(group_name) %>% summarise(pcb_cov = sum(width))
          
          ex_fo3 <- findOverlaps(IRanges::intersect(gr_blacklist, gr_test), exons_red)
          ex_fo_pint3 <- pintersect(IRanges::intersect(gr_blacklist, gr_test)[queryHits(ex_fo3)], exons_red[subjectHits(ex_fo3)])
          table_refseq_bl <- as.data.table(ex_fo_pint3) %>% group_by(group_name) %>% summarise(pcb_bl = sum(width))
          
          table_refseq <- left_join(
            left_join(table_refseq_cov, table_refseq_bl),
            table_refseq_covbl) %>% mutate(pcb_tot = sum(width(exons_red)), pcb_ncov = pcb_tot - pcb_covt, pcb_covp = (pcb_cov / pcb_tot)*100, pcb_covtp = (pcb_covt / pcb_tot)*100)
          colnames(table_refseq)[1] <- "gene"
          table_refseq <- table_refseq[,c("gene","pcb_ncov","pcb_covt","pcb_bl","pcb_cov","pcb_tot","pcb_covp","pcb_covtp")]
          
          # CLINVAR TABLE
          incProgress(0.1, detail = "Create ClinVar table")
          
          # create clinvar ranges and overlap with test ranges
          sqldb_clinvar <- dbConnect(SQLite(), dbname=clv_path)
          
          clvlabs <- unique(dbGetQuery(sqldb_clinvar, paste0('select clnsig from clinvar'))[,1])
          clvlabs <- clvlabs[str_starts(clvlabs, "Pathogenic") | str_starts(clvlabs, "Likely_pathogenic")]
          
          gr_clinvar <- GRanges(dbGetQuery(sqldb_clinvar, paste0('select CHROM from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                IRanges(dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                        dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1]+
                                          (nchar(dbGetQuery(sqldb_clinvar, paste0('select REF from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1])-1)))
          
          # total muts in targeted genes
          clv_gene_id <- dbGetQuery(sqldb_clinvar, paste0('select gene, ID from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))
          table_clv_tot <- as.data.frame(table(clv_gene_id$gene))
          
          # covered variants (excluding blacklisted)
          clv_fo_bl <- findOverlaps(gr_clinvar, gr_test_bl, type = "within")@from
          clv_cov_bl <- distinct(clv_gene_id[clv_fo_bl,])
          table_clv_covbl <- as.data.frame(sort(table(clv_cov_bl$gene)))
          
          # total and specific blacklisted variants
          clv_fo <- findOverlaps(gr_clinvar, gr_test, type = "within")@from
          clv_cov <- distinct(clv_gene_id[clv_fo,])
          table_clv_covt <- as.data.frame(table(clv_cov$gene))
          clv_bl <- clv_cov[!(clv_cov$ID %in% clv_cov_bl$ID),]
          #table_clv_bl <- as.data.frame(table(clv_bl$gene))
          if(!is.null(panelInput$mask)) {
            table_clv_bl <- as.data.frame(clv_bl$gene)
          } else {
            table_clv_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered variants (excluding blacklisted)
          clv_ncov <- clv_gene_id[!(clv_gene_id$ID %in% clv_cov$ID),]
          if (nrow(clv_ncov) > 0) {
            table_clv_ncov <- as.data.frame(table(clv_ncov$gene))
          } else {
            table_clv_ncov <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # join tables
          table_clv <- left_join(
            left_join(
              left_join(
                left_join(
                  full_join(data.frame("Var1" = gl_exp$refseq), table_clv_ncov, by = "Var1"),
                  table_clv_covt, by = "Var1"),
                table_clv_bl, by = "Var1"),
              table_clv_covbl, by = "Var1"),
            table_clv_tot, by = "Var1")
          colnames(table_clv) <- c("gene","clv_ncov","clv_covt","clv_bl","clv_cov","clv_tot")
          table_clv$clv_covp <- table_clv$clv_cov/table_clv$clv_tot
          table_clv$clv_covtp <- table_clv$clv_covt/table_clv$clv_tot
          table_clv[is.na(table_clv)] <- 0
          
          # Cosmic ------------------------------------------------------------------
          incProgress(0.3, detail = "Create COSMIC table")
          # make ranges, find overlap (blacklisted)
          
          sqldb_cosmic <- dbConnect(SQLite(), dbname=cmc_path)
          gr_cmc <-  GRanges(dbGetQuery(sqldb_cosmic, paste0('select chr from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                             IRanges(dbGetQuery(sqldb_cosmic, paste0('select start from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                     dbGetQuery(sqldb_cosmic, paste0('select end from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1]))
          
          
          # total muts in targeted genes
          cmc_gene_id <- dbGetQuery(sqldb_cosmic, paste0('select GENE_NAME, GENOMIC_MUTATION_ID, COSMIC_SAMPLE_POSRATE from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))
          table_cmc_tot <- as.data.frame(table(cmc_gene_id$GENE_NAME))
          
          # covered variants (excluding blacklisted)
          cmc_fo_bl <- findOverlaps(gr_cmc, gr_test_bl, type = "within")@from
          cmc_cov_bl <- distinct(cmc_gene_id[cmc_fo_bl,])
          table_cmc_covbl <- as.data.frame(table(cmc_cov_bl$GENE_NAME))
          
          # total and explicitly blacklisted variants
          cmc_fo <- findOverlaps(gr_cmc, gr_test, type = "within")@from
          cmc_cov <- distinct(cmc_gene_id[cmc_fo,])
          table_cmc_covt <- as.data.frame(table(cmc_cov$GENE_NAME))
          cmc_bl <- cmc_cov[!(cmc_cov$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          #table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
          if(!is.null(panelInput$mask)) {
            table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
          } else {
            table_cmc_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered variants (excluding blacklisted)
          cmc_ncov <- cmc_gene_id[!(cmc_gene_id$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
          if (nrow(cmc_ncov) > 0) {
            table_cmc_ncov <- as.data.frame(table(cmc_ncov$GENE_NAME))
          } else {
            table_cmc_ncov <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered including blacklisted
          cmc_ncovbl <- cmc_gene_id[!(cmc_gene_id$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          
          # rate of non-covered mutations - only in panel target genes
          cmc_ncov_posTestRate <- sum(cmc_ncov$COSMIC_SAMPLE_POSRATE)
          cmc_ncovbl_posTestRate <- sum(cmc_ncovbl$COSMIC_SAMPLE_POSRATE)
          
          cmc_id_posrate <- dbGetQuery(sqldb_cosmic, paste0('select GENOMIC_MUTATION_ID, COSMIC_SAMPLE_POSRATE from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other")'))
          
          # all non-covered mutations (including non-target genes)
          cmc_ncov_total <- cmc_id_posrate[!(cmc_id_posrate$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
          cmc_ncovbl_total <- cmc_id_posrate[!(cmc_id_posrate$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          
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
          panel <- list(panelName = input$pName, 
                        panelFullName = input$pfName,
                        panelTable = as_tibble(full_join(full_join(table_refseq, table_cmc, by = "gene"), table_clv, by = "gene")),
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
          
          dbxn$panelNames <<- c(dbxn$panelNames, input$pName)
          dbxn$geneNames <<- all_genes
          reset("pName")
          reset("pfName")
          reset("file1")
          reset("pRows")
          reset("pCols")
          reset("blackl")
          panelInput$data <- NULL
          panelInput$mask <- NULL
          incProgress(0.1, detail = "Cleaning up")
          gc()
          shinyalert(title = "complete")
        })
      }
    })
    
    
    # UPDATE ALL --------------------------------------------------------------
    
    observeEvent(input$updateb, {
      
      prepRefSeq(updateDb = T)
      loadRefSeq(force = T)
      prepClinvar(updateDb = updateCheckClv())
      prepCosmic()
      
      withProgress(message = "Updating", {
        for (j in 1:length(dbx)) {
          incProgress(1/length(dbx), detail = dbx[[j]][["panelName"]])
          gc()
          
          # PROCESS BED FILE
          #reorder file columns
          trget <- dbx[[j]][["panelBed_input"]]
          
          # make Granges of panel
          gr_test <- reduce(GRanges(trget$V1, IRanges(trget$V2, trget$V3)))
          
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
          gl_exp <- data.frame("refseq" = unique(names(unlist(ex_by_ge)[findOverlaps(gr_test, unlist(ex_by_ge))@to ])))
          
          # identify unique exons and their coverage
          exons <- ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq]
          ex_fo <- findOverlaps(gr_test, exons)
          ex_fo_bl <- findOverlaps(gr_test_bl, exons)
          ex_by_ge_all_df <- subset(
            left_join(as.data.table(exons), left_join(
              as.data.table(pintersect(gr_test[queryHits(ex_fo)], exons[subjectHits(ex_fo)])) %>% 
                filter(hit == T) %>%
                group_by(exon_id) %>%
                summarise(covt_width = sum(width)),
              as.data.table(pintersect(gr_test_bl[queryHits(ex_fo_bl)], exons[subjectHits(ex_fo_bl)])) %>% 
                filter(hit == T) %>%
                group_by(exon_id) %>%
                summarise(cov_width = sum(width))),
                      by = "exon_id"),
            select = c(group_name, seqnames, start, end, strand, width, covt_width, cov_width)) %>%
            distinct() %>%
            replace(is.na(.), 0)
          
          # reduce exon ranges
          exons_red <- reduce(ex_by_ge[names(ex_by_ge) %in% gl_exp$refseq])
          
          ex_fo1 <- findOverlaps(gr_test, exons_red)
          ex_fo_pint1 <- pintersect(gr_test[queryHits(ex_fo1)], exons_red[subjectHits(ex_fo1)])
          table_refseq_cov <- as.data.table(ex_fo_pint1) %>% group_by(group_name) %>% summarise(pcb_covt = sum(width))
          
          ex_fo2 <- findOverlaps(gr_test_bl, exons_red)
          ex_fo_pint2 <- pintersect(gr_test_bl[queryHits(ex_fo2)], exons_red[subjectHits(ex_fo2)])
          table_refseq_covbl <- as.data.table(ex_fo_pint2) %>% group_by(group_name) %>% summarise(pcb_cov = sum(width))
          
          ex_fo3 <- findOverlaps(IRanges::intersect(gr_blacklist, gr_test), exons_red)
          ex_fo_pint3 <- pintersect(IRanges::intersect(gr_blacklist, gr_test)[queryHits(ex_fo3)], exons_red[subjectHits(ex_fo3)])
          table_refseq_bl <- as.data.table(ex_fo_pint3) %>% group_by(group_name) %>% summarise(pcb_bl = sum(width))
          
          table_refseq <- left_join(
            left_join(table_refseq_cov, table_refseq_bl),
            table_refseq_covbl) %>% mutate(pcb_tot = sum(width(exons_red)), pcb_ncov = pcb_tot - pcb_covt, pcb_covp = (pcb_cov / pcb_tot)*100, pcb_covtp = (pcb_covt / pcb_tot)*100)
          colnames(table_refseq)[1] <- "gene"
          table_refseq <- table_refseq[,c("gene","pcb_ncov","pcb_covt","pcb_bl","pcb_cov","pcb_tot","pcb_covp","pcb_covtp")]
          
          # CLINVAR TABLE
          
          # create clinvar ranges and overlap with test ranges
          sqldb_clinvar <- dbConnect(SQLite(), dbname=clv_path)
          
          clvlabs <- unique(dbGetQuery(sqldb_clinvar, paste0('select clnsig from clinvar'))[,1])
          clvlabs <- clvlabs[str_starts(clvlabs, "Pathogenic") | str_starts(clvlabs, "Likely_pathogenic")]
          
          gr_clinvar <- GRanges(dbGetQuery(sqldb_clinvar, paste0('select CHROM from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                IRanges(dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                        dbGetQuery(sqldb_clinvar, paste0('select POS from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1]+
                                          (nchar(dbGetQuery(sqldb_clinvar, paste0('select REF from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1])-1)))
          
          # total muts in targeted genes
          clv_gene_id <- dbGetQuery(sqldb_clinvar, paste0('select gene, ID from clinvar where (clnsig in (', paste(shQuote(clvlabs, type = "cmd"), collapse = ", "),') and gene in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))
          table_clv_tot <- as.data.frame(table(clv_gene_id$gene))
          
          # covered variants (excluding blacklisted)
          clv_fo_bl <- findOverlaps(gr_clinvar, gr_test_bl, type = "within")@from
          clv_cov_bl <- distinct(clv_gene_id[clv_fo_bl,])
          table_clv_covbl <- as.data.frame(sort(table(clv_cov_bl$gene)))
          
          # total and specific blacklisted variants
          clv_fo <- findOverlaps(gr_clinvar, gr_test, type = "within")@from
          clv_cov <- distinct(clv_gene_id[clv_fo,])
          table_clv_covt <- as.data.frame(table(clv_cov$gene))
          clv_bl <- clv_cov[!(clv_cov$ID %in% clv_cov_bl$ID),]
          #table_clv_bl <- as.data.frame(table(clv_bl$gene))
          if(!is.null(panelInput$mask)) {
            table_clv_bl <- as.data.frame(clv_bl$gene)
          } else {
            table_clv_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered variants (excluding blacklisted)
          clv_ncov <- clv_gene_id[!(clv_gene_id$ID %in% clv_cov$ID),]
          if (nrow(clv_ncov) > 0) {
            table_clv_ncov <- as.data.frame(table(clv_ncov$gene))
          } else {
            table_clv_ncov <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # join tables
          table_clv <- left_join(
            left_join(
              left_join(
                left_join(
                  full_join(data.frame("Var1" = gl_exp$refseq), table_clv_ncov, by = "Var1"),
                  table_clv_covt, by = "Var1"),
                table_clv_bl, by = "Var1"),
              table_clv_covbl, by = "Var1"),
            table_clv_tot, by = "Var1")
          colnames(table_clv) <- c("gene","clv_ncov","clv_covt","clv_bl","clv_cov","clv_tot")
          table_clv$clv_covp <- table_clv$clv_cov/table_clv$clv_tot
          table_clv$clv_covtp <- table_clv$clv_covt/table_clv$clv_tot
          table_clv[is.na(table_clv)] <- 0
          
          # Cosmic ------------------------------------------------------------------
          # make ranges, find overlap (blacklisted)
          
          sqldb_cosmic <- dbConnect(SQLite(), dbname=cmc_path)
          gr_cmc <-  GRanges(dbGetQuery(sqldb_cosmic, paste0('select chr from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                             IRanges(dbGetQuery(sqldb_cosmic, paste0('select start from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1],
                                     dbGetQuery(sqldb_cosmic, paste0('select end from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))[,1]))
          
          
          # total muts in targeted genes
          cmc_gene_id <- dbGetQuery(sqldb_cosmic, paste0('select GENE_NAME, GENOMIC_MUTATION_ID, COSMIC_SAMPLE_POSRATE from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other" and GENE_NAME in (', paste(shQuote(gl_exp$refseq, type = "cmd"), collapse = ", "),'))'))
          table_cmc_tot <- as.data.frame(table(cmc_gene_id$GENE_NAME))
          
          # covered variants (excluding blacklisted)
          cmc_fo_bl <- findOverlaps(gr_cmc, gr_test_bl, type = "within")@from
          cmc_cov_bl <- distinct(cmc_gene_id[cmc_fo_bl,])
          table_cmc_covbl <- as.data.frame(table(cmc_cov_bl$GENE_NAME))
          
          # total and explicitly blacklisted variants
          cmc_fo <- findOverlaps(gr_cmc, gr_test, type = "within")@from
          cmc_cov <- distinct(cmc_gene_id[cmc_fo,])
          table_cmc_covt <- as.data.frame(table(cmc_cov$GENE_NAME))
          cmc_bl <- cmc_cov[!(cmc_cov$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          #table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
          if(!is.null(panelInput$mask)) {
            table_cmc_bl <- as.data.frame(table(cmc_bl$GENE_NAME))
          } else {
            table_cmc_bl <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered variants (excluding blacklisted)
          cmc_ncov <- cmc_gene_id[!(cmc_gene_id$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
          if (nrow(cmc_ncov) > 0) {
            table_cmc_ncov <- as.data.frame(table(cmc_ncov$GENE_NAME))
          } else {
            table_cmc_ncov <- data.frame("Var1" = gl_exp$refseq, "Freq" = 0)
          }
          
          # NON covered including blacklisted
          cmc_ncovbl <- cmc_gene_id[!(cmc_gene_id$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          
          # rate of non-covered mutations - only in panel target genes
          cmc_ncov_posTestRate <- sum(cmc_ncov$COSMIC_SAMPLE_POSRATE)
          cmc_ncovbl_posTestRate <- sum(cmc_ncovbl$COSMIC_SAMPLE_POSRATE)
          
          cmc_id_posrate <- dbGetQuery(sqldb_cosmic, paste0('select GENOMIC_MUTATION_ID, COSMIC_SAMPLE_POSRATE from cosmic where (MUTATION_SIGNIFICANCE_TIER != "Other")'))
          
          # all non-covered mutations (including non-target genes)
          cmc_ncov_total <- cmc_id_posrate[!(cmc_id_posrate$GENOMIC_MUTATION_ID %in% cmc_cov$GENOMIC_MUTATION_ID),]
          cmc_ncovbl_total <- cmc_id_posrate[!(cmc_id_posrate$GENOMIC_MUTATION_ID %in% cmc_cov_bl$GENOMIC_MUTATION_ID),]
          
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
                        panelTable = as_tibble(full_join(full_join(table_refseq, table_cmc, by = "gene"), table_clv, by = "gene")),
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
