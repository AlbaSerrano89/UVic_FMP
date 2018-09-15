function(input, output, session)
{
  # Loading the functions we will need
  source("./Functions.R")
  
  # Reading the microarray data, with, at least, one column of Entrez IDs, another one with the pvalues and the last one with the log Fold Change values.
  csv <- reactive({
    inFile <- input$microarray
    
    csv <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
    csv
  })
  
  # Selecting the differentially expressed genes. "pval_init" and "logfc" are values written by the user.
  deg <- reactive({
    csv <- csv()
    
    pval_col <- which(names(csv) == input$cols_pval)
    logfc_col <- which(names(csv) == input$cols_logfc)

    deg <- csv[csv[, pval_col] < input$pval_init &
                 (csv[, logfc_col] > input$logfc |
                    csv[, logfc_col] < -input$logfc), ]
    deg
  })
  
  # The initial data could have more than one gene per transcript.
  #     If the user decides to remove these cases, we get "rem" from the input.
  #     If the user decides to use all the genes, we get "use" from the input.
  rem_use <- reactive({
    deg <- deg()
    entrez_col <- which(names(deg) == input$cols_entrez)
    
    if(input$numb_genes == "rem")
    {
      rem <- which(!is.na(as.numeric(as.character(deg[, entrez_col]))))
      rem_use <- deg[rem, ]
    }
    else rem_use <- repetitions(deg, entrez_col, " /// ")
    
    rem_use
  })
  
  # Reading the database we have created with Uniprot about the specie the user wants to analyse.
  res <- reactive({
    if(input$specie == 9606)
    {
      file <- "data_HSapiens.csv"
      variables <- input$human_var
    }
    else if(input$specie == 10090)
    {
      file <- "data_MMusculus.csv"
      variables <- input$mouse_var
    }
    else if(input$specie == 3702)
    {
      file <- "data_AThaliana.csv"
      variables <- input$cress_var
    }
    else if(input$specie == 7227)
    {
      file <- "data_DMelanogaster.csv"
      variables <- input$fly_var
      library("drosophila2.db")
    }
    else if(input$specie == 6239)
      {
      file <- "data_CElegans.csv"
      variables <- input$worm_var
    }
    else if(input$specie == 559292)
    {
      file <- "data_SCerevisiae.csv"
      variables <- input$yeast_var
    }
    
    res <- read.csv(file, stringsAsFactors = FALSE)
    res <- res[!is.na(res$UNIPROTKB), c(1, 2, as.numeric(variables))]
    res
  })
  
  # Creating a big dataframe, by joining the differentially expressed genes got from the initial data with the data of the Uniprot.
  # We only want one protein per row, so we join the values with "; " if there are more than one value in some field.
  big_dataframe <- reactive({
    rem_use <- rem_use()
    res <- res()
    entrez_col <- which(names(rem_use) == input$cols_entrez)
    
    de_genes <- rem_use[, entrez_col]
    de_prots_i <- NULL
    for(i in 1:length(de_genes))
    {
      a <- grep(de_genes[i], res$ENTREZ_GENE)
      if(length(a) > 0) de_prots_i <- c(de_prots_i, a)
    }
    
    deg_prots <- res[unique(de_prots_i), ]
    
    fill <- as.data.frame(matrix(rep("NA", nrow(deg_prots) * ncol(rem_use)), nrow = nrow(deg_prots)), stringsAsFactors = FALSE)
    for(i in 1:nrow(deg_prots))
    {
      a <- grep(deg_prots$ENTREZ_GENE[i], rem_use[, entrez_col])
      if(length(a) == 1) fill[i, ] <- rem_use[a, ]
      else
      {
        b <- NULL
        for(j in 1:ncol(rem_use))
        {
          f <- paste(unique(rem_use[a, j]), collapse = "; ")
          b <- c(b, f)
        }
        fill[i, ] <- b
      }
    }
    colnames(fill) <- colnames(rem_use)
    
    big_dataframe <- cbind(deg_prots, fill)
    big_dataframe
  })
  
  # Reading the predictions obtained by TPpred2.0
  outputs_TPpred <- reactive({
    if(input$specie == 9606) file <- "TPpred_HSapiens.txt"
    else if(input$specie == 10090) file <- "TPpred_MMusculus.txt"
    else if(input$specie == 3702) file <- "TPpred_AThaliana.txt"
    else if(input$specie == 7227) file <- "TPpred_DMelanogaster.txt"
    else if(input$specie == 6239) file <- "TPpred_CElegans.txt"
    else if(input$specie == 559292) file <- "TPpred_SCerevisiae.txt"
    
    outputs_TPpred <- total_percs(file)
    outputs_TPpred
  })
  
  # Joining the predictions of the TPpred2.0 with the big dataframe created before.
  all_data <- reactive({
    outputs_TPpred <- outputs_TPpred()
    big_dataframe <- big_dataframe()
    
    colnames(outputs_TPpred)[1] <- colnames(big_dataframe)[1]
    data_prot <- merge(big_dataframe, outputs_TPpred, by = intersect(colnames(big_dataframe), colnames(outputs_TPpred)))
    data_prot <- data_prot[-which(data_prot$Probe.Set.ID == ""), ]
    
    logfc_dataprot_i <- which(names(data_prot) == input$cols_logfc)
    
    up_down_reg <- NULL
    for(i in 1:nrow(data_prot))
    {
      if(length(grep("; ", data_prot[i, logfc_dataprot_i])) == 0)
      {
        if(data_prot[i, logfc_dataprot_i] < 0) up_down_reg[i] <- "down"
        else up_down_reg[i] <- "up"
      }
      else
      {
        logs <- unlist(strsplit(data_prot[i, logfc_dataprot_i], split = "; "))
        a <- NULL
        for(j in 1:length(logs))
        {
          if(as.numeric(logs)[j] < 0) a[j] <- "down"
          else a[j] <- "up"
        }
        if(sum(duplicated(a)) == (length(logs) - 1)) up_down_reg[i] <- a[1]
        else up_down_reg[i] <- "up_down"
      }
    }
    
    all_data1 <- data.frame(up_down_reg, data_prot)
    ncols_mega <- length(all_data1)
    pval_col <- which(names(all_data1) == input$cols_pval)
    all_data2 <- all_data1[, c(2, ncols_mega - 1, ncols_mega, 3, 1, 4:(ncols_mega - 2))]
    
    pval_all_data_i <- which(names(all_data2) == input$cols_pval)
    new_pvals <-NULL
    for(i in 1:nrow(all_data2))
    {
      if(length(grep(";", all_data2$pval[i])) == 0) new_pvals[i] <- as.numeric(all_data2[i, pval_all_data_i])
      else new_pvals[i] <- mean(as.numeric(unlist(strsplit(all_data2[i, pval_all_data_i], split = "; "))))
    }
    
    all_data <- all_data2[order(new_pvals), ]
    all_data
  })
  
  # The same as before, but just taking those genes that are "up" or "down", but not both.
  mega_dataframe <- reactive({
    all_data <- all_data()
    mega_dataframe <- all_data[all_data$up_down_reg != "up_down", ]
    mega_dataframe
  })
  
  # The same as before, but just taking those genes that are not "up" nor "down".
  up_down_dataframe <- reactive({
    all_data <- all_data()
    up_down_dataframe <- all_data[all_data$up_down_reg == "up_down", ]
    up_down_dataframe
  })
  
  # Preparing the data to do the GO analysis: only entrez, mito_status, up_down_reg & GO_ID columns; one gene per row.
  finalGenesMitoDE <- reactive({
    mega_dataframe <- mega_dataframe()
    entrez_col_mega <- which(names(mega_dataframe) == input$cols_entrez)
    
    all_emde <- mega_dataframe[, c(entrez_col_mega, 2, 5, 10)]
    all_emde$up_down_reg <- as.character(all_emde$up_down_reg)
    
    entrez_col_emde <- which(names(all_emde) == "Entrez.Gene")
    all_emde2 <- all_emde[order(as.numeric(all_emde[, entrez_col_emde])), ]
    all_emde2 <- all_emde2[!duplicated(all_emde2), ]
    entrez_col_emde2 <- which(names(all_emde2) == "Entrez.Gene")
    
    mitos <- NULL
    for(i in 1:nrow(all_emde2))
    {
      all_i <- all_emde2[which(all_emde2[, entrez_col_emde2] == all_emde2[i, entrez_col_emde2]), ]
      
      if(nrow(all_i) == 1) mitos[i] <- as.character(all_i[1, 2])
      else
      {
        if(sum(duplicated(all_i[, 2])) == (nrow(all_i) - 1)) mitos[i] <- as.character(all_i[1, 2])
        else mitos[i] <- "undetermined"
      }
    }
    all_emde2$mito_status <- mitos
    
    finalGenesMitoDE <- all_emde2[!duplicated(all_emde2[, entrez_col_emde2]), ]
    finalGenesMitoDE
  })
  
  # Downloading the GO names of the differentially expressed genes, from Uniprot.
  ent_go <- reactive({
    finalGenesMitoDE <- finalGenesMitoDE()
    
    entrez_col_fGenesMDE <- which(names(finalGenesMitoDE) == input$cols_entrez)
    
    library(UniProt.ws)
    unip <- UniProt.ws(taxId = as.numeric(input$specie))
    ent_go <- select(unip,
                     keys = as.numeric(finalGenesMitoDE[, entrez_col_fGenesMDE]),
                     columns = c("ENTREZ_GENE", "GO-ID", "GO"),
                     keytype = "ENTREZ_GENE")
    colnames(ent_go)[2] <- "GO.ID"
    ent_go
  })
  
  # Starting to create a dataframe with GO names as rows (one GO per row), and the rest of information from "finalGenesMitoDE" (1/4).
  newGOs <- reactive({
    ent_go <- ent_go()
    finalGenesMitoDE <- finalGenesMitoDE()
    
    gos <- strsplit(finalGenesMitoDE$GO.ID, split = "; ")
    
    GOnames <- NULL
    for(i in 1:nrow(finalGenesMitoDE))
    {
      GOi <- NULL
      for(j in 1:length(gos[[i]]))
      {
        quins <- grep(gos[[i]][[j]], ent_go$GO)
        gous <- unlist(strsplit(ent_go$GO[quins[1]], split = "; "))
        GOname <- gous[grep(gos[[i]][[j]], gous)]
        GOi[j] <- GOname
      }
      GOnames[i] <- paste(GOi, collapse = "; ")
    }
    finalGenesMitoDE$GO <- GOnames
    newGOs <- finalGenesMitoDE
    newGOs <- newGOs[-which(is.na(newGOs$GO.ID)), ]
    newGOs
  })
  
  # Table: number of genes per GO of the new data.
  FinalNewGOs <- reactive({
    newGOs <- newGOs()
    
    newGOs_res <- strsplit(newGOs$GO.ID, split = "; ")
    entrez_col_newGOs <- which(names(newGOs) == input$cols_entrez)
    names(newGOs_res) <- newGOs[, entrez_col_newGOs]
    
    FinalNewGOs <- sort(table(unlist(newGOs_res)), decreasing = TRUE)
    FinalNewGOs
  })
  
  # Proportion of genes per GO of the new data.
  numberNewGO <- reactive({
    newGOs <- newGOs()
    FinalNewGOs <- FinalNewGOs()
    
    entrez_col_newGOs <- which(names(newGOs) == input$cols_entrez)
    
    newGenes <- length(unique(newGOs[, entrez_col_newGOs]))
    numberNewGO <- FinalNewGOs / newGenes * 100
    numberNewGO
  })
  
  # Continuing creating a dataframe with GO names as rows (one GO per row), and the rest of information from "finalGenesMitoDE" (2/4).
  dataset_new <- reactive({
    FinalNewGOs <- FinalNewGOs()
    newGOs <- newGOs()
    
    entrez_col_newGOs <- which(names(newGOs) == input$cols_entrez)
    
    dataset_new <- data.frame(character(length(FinalNewGOs)),
                              character(length(FinalNewGOs)),
                              character(length(FinalNewGOs)),
                              stringsAsFactors = FALSE)
    colnames(dataset_new) <- c("GO_id", "TOTAL", "AllGenes")
    for(i in 1:length(FinalNewGOs))
    {
      quins <- grep(names(FinalNewGOs[i]), newGOs$GO.ID)
      new_row <- paste(newGOs[quins, entrez_col_newGOs], collapse = "; ")
      dataset_new[i, ] <- c(names(FinalNewGOs)[i], as.numeric(FinalNewGOs[i]), new_row)
    }
    
    namesGO <- unique(unlist(strsplit(newGOs$GO, split = "; ")))
    colGOnames <- NULL
    for(i in 1:nrow(dataset_new))
    {
      go <- grep(dataset_new$GO_id[i], namesGO)
      colGOnames[i] <- namesGO[go[1]]
    }
    dataset_new$GO <- colGOnames
    
    dataset_new
  })
  
  # Reading the dataset about the GO terms.
  totalGOs <- reactive({
    if(input$specie == 9606) file <- "resGO_HS.csv"
    else if(input$specie == 10090) file <- "resGO_MM.csv"
    else if(input$specie == 3702) file <- "resGO_AT.csv"
    else if(input$specie == 7227) file <- "resGO_DM.csv"
    else if(input$specie == 6239) file <- "resGO_CE.csv"
    else if(input$specie == 559292) file <- "resGO_SC.csv"
    
    totalGOs <- read.csv(file, stringsAsFactors = FALSE)
    totalGOs <- totalGOs[-which(totalGOs$GO.ID == ""), ]
    totalGOs
  })
  
  # Table: number of genes per GO of the whole genome.
  FinalAllGOs <- reactive({
    FinalNewGOs <- FinalNewGOs()
    totalGOs <- totalGOs()
    
    FinalAllGOs <- NULL
    for(i in 1:length(FinalNewGOs))
    {
      go_i <- totalGOs[grep(names(FinalNewGOs)[i], totalGOs$GO.ID), ]
      FinalAllGOs[i] <- nrow(go_i)
    }
    names(FinalAllGOs) <- names(FinalNewGOs)
    FinalAllGOs
  })
  
  # Proportion of genes per GO of the whole genome.
  numberAllGO <- reactive({
    totalGOs <- totalGOs()
    FinalAllGOs <- FinalAllGOs()
    
    totalGenes <- length(totalGOs$ENTREZ_GENE)
    numberAllGO <- FinalAllGOs / totalGenes * 100
    numberAllGO
  })
  
  # Performing one proportion test per GO term, comparing the DEG to the whole genome, and adding the pvalues to the dataframe (3/4).
  prop_test_GO <- reactive({
    dataset_new <- dataset_new()
    FinalNewGOs <- FinalNewGOs()
    FinalAllGOs <- FinalAllGOs()
    newGOs <- newGOs()
    totalGOs <- totalGOs()
    
    entrez_col_newGOs <- which(names(newGOs) == input$cols_entrez)
    
    newGenes <- length(unique(newGOs[, entrez_col_newGOs]))
    totalGenes <- length(totalGOs$ENTREZ_GENE)
    
    significance <- (100 - input$confidenceGO) / 100
    
    pval_testGO <- NULL
    significative <- NULL
    for(i in 1:nrow(dataset_new))
    {
      test <- prop.test(c(FinalNewGOs[i], FinalAllGOs[i]), c(newGenes, totalGenes), alt = "t")
      pval_testGO[i] <- test$p.value
      
      if(test$p.value >= significance) significative[i] <- "no"
      else
      {
        if(FinalNewGOs[i] / newGenes > FinalAllGOs[i] / totalGenes) significative[i] <- "enriched"
        else significative[i] <- "depleted"
      }
    }
    
    list(pval_testGO, significative)
  })
  
  # Finalizing the GO dataset (4/4).
  DatasetGO <- reactive({
    dataset_new <- dataset_new()
    prop_test_GO <- prop_test_GO()
    ent_go <- ent_go()
    newGOs <- newGOs()
    
    entrez_col_newGOs <- which(names(newGOs) == input$cols_entrez)
    
    dataset_new$PvalTestProp <- prop_test_GO[[1]]
    dataset_new$Significative <- prop_test_GO[[2]]
    
    genes <- strsplit(as.character(dataset_new$AllGenes), split = "; ")
    names(genes) <- dataset_new$GO_id
    
    new_dataset_GO <- data.frame(character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 character(length(genes)),
                                 stringsAsFactors = FALSE)
    colnames(new_dataset_GO) <- c("Mito_Down", "Mito_Up", "NoMito_Down", "NoMito_Up", "Genes-Mito_Down", "Genes-Mito_Up", "Genes-NoMito_Down", "Genes-NoMito_Up")
    for(i in 1:length(genes))
    {
      GOterm <- genes[[i]]
      red_df <- newGOs[newGOs[, entrez_col_newGOs] %in% GOterm, ]
      
      md_i <- which(red_df$mito_status == "mito" & red_df$up_down_reg == "down")
      mu_i <- which(red_df$mito_status == "mito" & red_df$up_down_reg == "up")
      nd_i <- which(red_df$mito_status == "no_mito" & red_df$up_down_reg == "down")
      nu_i <- which(red_df$mito_status == "no_mito" & red_df$up_down_reg == "up")
      
      md_genes <- paste(red_df[md_i, 1], collapse = "; ")
      mu_genes <- paste(red_df[mu_i, 1], collapse = "; ")
      nd_genes <- paste(red_df[nd_i, 1], collapse = "; ")
      nu_genes <- paste(red_df[nu_i, 1], collapse = "; ")
      
      new_dataset_GO[i, ] <- c(length(md_i), length(mu_i), length(nd_i), length(nu_i),
                               md_genes, mu_genes, nd_genes, nu_genes)
    }
    
    final_dataset <- cbind(dataset_new, new_dataset_GO)
    final_dataset <- final_dataset[, c(4:10, 2, 11:14, 1, 3)]
    DatasetGO <- final_dataset[order(final_dataset$PvalTestProp), ]
    DatasetGO
  })
  
  # Table: number of mito and no_mito proteins of the input data.
  nums_new_mito <- reactive({
    mega_dataframe <- mega_dataframe()
    nums_new_mito <- addmargins(table(mega_dataframe$mito_status))
    names(nums_new_mito) <- c("mito", "no_mito", "total")
    nums_new_mito
  })
  
  # Percentages of mito and no_mito proteins of the input data.
  perc_new_mito <- reactive({
    nums_new_mito <- nums_new_mito()
    perc_new_mito <- round(nums_new_mito[1:2] * 100 / nums_new_mito[3], 2)
    perc_new_mito
  })
  
  # Table: number of mito and no_mito proteins of the whole proteome.
  nums_total_mito <- reactive({
    outputs_TPpred <- outputs_TPpred()
    nums_total_mito <- addmargins(table(outputs_TPpred$mito_status))
    names(nums_total_mito) <- c("mito", "no_mito", "total")
    nums_total_mito
  })
  
  # Percentages of mito and no_mito proteins of the whole proteome.
  perc_total_mito <- reactive({
    nums_total_mito <- nums_total_mito()
    perc_total_mito <- round(nums_total_mito[1:2] * 100 / nums_total_mito[3], 2)
    perc_total_mito
  })
  
  # Table: number of mito and no_mito proteins, crossed with number of up and down genes of the input data.
  tab_deg <- reactive({
    mega_dataframe <- mega_dataframe()
    table(mega_dataframe$mito_status, mega_dataframe$up_down_reg)
  })
  
  # Vector of the last table.
  nums_tab_deg <- reactive({
    tab_deg <- tab_deg()
    
    nums_tab_deg <- c(tab_deg)
    nums_tab_deg
  })
  
  # Percentages of mito and no_mito proteins, crossed with number of up and down genes of the input data.
  perc_tab_deg <- reactive({
    nums_tab_deg <- nums_tab_deg()

    ntotal <- sum(nums_tab_deg)
    perc_tab_deg <- round(nums_tab_deg * 100 / ntotal, 2)
    perc_tab_deg
  })
  
  # Performing the mito proportion test, comparing the mito results of the differentially expressed proteins to the whole proteome ones.
  # Giving the conclusions of the test.
  conclusions <- reactive({
    nums_total_mito <- nums_total_mito()
    nums_new_mito <- nums_new_mito()
    significance <- (100 - input$confidence) / 100
    
    test_t <- prop.test(c(nums_new_mito[1], nums_total_mito[1]), c(nums_new_mito[3], nums_total_mito[3]), alt = "t")
    
    if(test_t$p.value >= significance)
    {
      pval_test <- paste("The p-value of the proportion test is", round(test_t$p.value, 4))
      conc <- paste0("With a confidence of ", input$confidence, "%,\nthere are NOT enough evidences to conclude that there are significant differences between your dataset and the whole proteome.")
    }
    else
    {
      if(nums_new_mito[1] / nums_new_mito[3] > nums_total_mito[1] / nums_total_mito[3])
      {
        test_g <- prop.test(c(nums_new_mito[1], nums_total_mito[1]), c(nums_new_mito[3], nums_total_mito[3]), alt = "g")
        pval_test <- paste("The p-value of the proportion test is", round(test_g$p.value, 4))
        conc <- paste0("With a confidence of ", input$confidence, "%,\nthere are enough evidences to conclude that your dataset is significantly enriched in mitochondrial sequences, regarding the whole proteome.")
      }
      else
      {
        test_l <- prop.test(c(nums_new_mito[1], nums_total_mito[1]), c(nums_new_mito[3], nums_total_mito[3]), alt = "l")
        pval_test <- paste("The p-value of the proportion test is", round(test_l$p.value, 4))
        conc <- paste0("With a confidence of ", input$confidence, "%,\nthere are enough evidences to conclude that your dataset is significantly depleted in mitochondrial sequences, regarding the whole proteome.")
      }
    }
    
    list(pval_test, conc)
  })
  
  # Drawing the Volcano plot of the DEG.
  volc <- reactive({
    csv <- csv()
    
    pval_col <- which(names(csv) == input$cols_pval)
    logfc_col <- which(names(csv) == input$cols_logfc)
    
    lim_logfc <- input$logfc
    lim_pval <- input$pval_init
    
    x_lim <- max(abs(min(csv[, logfc_col]) - 0.5), max(csv[, logfc_col]) + 0.5)
    y_lim <- max(-log10(csv[, pval_col])) + 0.5
    
    sel_up <- csv[csv[, pval_col] < lim_pval & csv[, logfc_col] > lim_logfc, ]
    sel_do <- csv[csv[, pval_col] < lim_pval & csv[, logfc_col] < -lim_logfc, ]
    not_sel1 <- csv[csv[, pval_col] >= lim_pval & (csv[, logfc_col] > -lim_logfc & csv[, logfc_col] < lim_logfc), ]
    not_sel2 <- csv[csv[, pval_col] < lim_pval & (csv[, logfc_col] > -lim_logfc & csv[, logfc_col] < lim_logfc), ]
    
    plot(not_sel1[, logfc_col], -log10(not_sel1[, pval_col]),
         pch = 20, cex = 0.8,
         main = "Volcano plot",
         xlab = "log2(FC)", ylab = "-log10(pvalue)",
         xlim = c(-x_lim, x_lim), ylim = c(0, y_lim))
    points(not_sel2[, logfc_col], -log10(not_sel2[, pval_col]),
           pch = 20, cex = 0.8)
    points(sel_up[, logfc_col], -log10(sel_up[, pval_col]),
           pch = 20, cex = 0.8,
           col = "red")
    points(sel_do[, logfc_col], -log10(sel_do[, pval_col]),
           pch = 20, cex = 0.8,
           col = "blue")
    abline(h = -log10(lim_pval), col = "gray")
    abline(v = lim_logfc, col = "gray")
    abline(v = -lim_logfc, col = "gray")
  })
  
  # Drawing the Volcano plot of the DEG. (download)
  volc2 <- function(){
    csv <- csv()
    
    pval_col <- which(names(csv) == input$cols_pval)
    logfc_col <- which(names(csv) == input$cols_logfc)
    
    lim_logfc <- input$logfc
    lim_pval <- input$pval_init
    
    x_lim <- max(abs(min(csv[, logfc_col]) - 0.5), max(csv[, logfc_col]) + 0.5)
    y_lim <- max(-log10(csv[, pval_col])) + 0.5
    
    sel_up <- csv[csv[, pval_col] < lim_pval & csv[, logfc_col] > lim_logfc, ]
    sel_do <- csv[csv[, pval_col] < lim_pval & csv[, logfc_col] < -lim_logfc, ]
    not_sel1 <- csv[csv[, pval_col] >= lim_pval & (csv[, logfc_col] > -lim_logfc & csv[, logfc_col] < lim_logfc), ]
    not_sel2 <- csv[csv[, pval_col] < lim_pval & (csv[, logfc_col] > -lim_logfc & csv[, logfc_col] < lim_logfc), ]
    
    plot(not_sel1[, logfc_col], -log10(not_sel1[, pval_col]),
         pch = 20, cex = 0.8,
         main = "Volcano plot",
         xlab = "log2(FC)", ylab = "-log10(pvalue)",
         xlim = c(-x_lim, x_lim), ylim = c(0, y_lim))
    points(not_sel2[, logfc_col], -log10(not_sel2[, pval_col]),
           pch = 20, cex = 0.8)
    points(sel_up[, logfc_col], -log10(sel_up[, pval_col]),
           pch = 20, cex = 0.8,
           col = "red")
    points(sel_do[, logfc_col], -log10(sel_do[, pval_col]),
           pch = 20, cex = 0.8,
           col = "blue")
    abline(h = -log10(lim_pval), col = "gray")
    abline(v = lim_logfc, col = "gray")
    abline(v = -lim_logfc, col = "gray")
  }
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins, where the total are the differentially expressed ones.
  pieNew <- reactive({
    perc_new_mito <- perc_new_mito()
    nums_new_mito <- nums_new_mito()
    
    percs_lbls <- paste0(perc_new_mito, "%")
    
    lab1 <- paste0("Predicted Mitochondrial Location: ", nums_new_mito[1])
    lab2 <- paste0("Predicted Ex-Mitochondrial Location: ", nums_new_mito[2])
    lbls <- c(lab1, lab2)
    
    main_pie <- paste("YOUR DATA\nTotal proteins:", nums_new_mito[3])
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(perc_new_mito, radius = 0.8, labels = percs_lbls, col = c("red", "cyan"))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "cyan"))
  })
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins, where the total are the differentially expressed ones. (download)
  pieNew2 <- function(){
    perc_new_mito <- perc_new_mito()
    nums_new_mito <- nums_new_mito()
    
    percs_lbls <- paste0(perc_new_mito, "%")
    
    lab1 <- paste0("Predicted Mitochondrial Location: ", nums_new_mito[1])
    lab2 <- paste0("Predicted Ex-Mitochondrial Location: ", nums_new_mito[2])
    lbls <- c(lab1, lab2)
    
    main_pie <- paste("YOUR DATA\nTotal proteins:", nums_new_mito[3])
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(perc_new_mito, radius = 0.8, labels = percs_lbls, col = c("red", "cyan"))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "cyan"))
  }
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins, where the total are the whole proteome ones.
  pieTotal <- reactive({
    perc_total_mito <- perc_total_mito()
    nums_total_mito <- nums_total_mito()
    
    percs_lbls <- paste0(perc_total_mito, "%")
    
    lab1 <- paste0("Predicted Mitochondrial Location: ", nums_total_mito[1])
    lab2 <- paste0("Predicted Ex-Mitochondrial Location: ", nums_total_mito[2])
    lbls <- c(lab1, lab2)
    
    main_pie <- paste("THE WHOLE PROTEOME\nTotal proteins:", nums_total_mito[3])
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(perc_total_mito, radius = 0.8, labels = percs_lbls, col = c("red", "cyan"))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "cyan"))
  })
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins, where the total are the whole proteome ones. (download)
  pieTotal2 <- function(){
    perc_total_mito <- perc_total_mito()
    nums_total_mito <- nums_total_mito()
    
    percs_lbls <- paste0(perc_total_mito, "%")
    
    lab1 <- paste0("Predicted Mitochondrial Location: ", nums_total_mito[1])
    lab2 <- paste0("Predicted Ex-Mitochondrial Location: ", nums_total_mito[2])
    lbls <- c(lab1, lab2)
    
    main_pie <- paste("THE WHOLE PROTEOME\nTotal proteins:", nums_total_mito[3])
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(perc_total_mito, radius = 0.8, labels = percs_lbls, col = c("red", "cyan"))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "cyan"))
  }
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins AND the proportions of up and down regulated genes, where the total are differentially expressed ones.
  pieDeg <- reactive({
    nums_tab_deg <- nums_tab_deg()
    perc_tab_deg <- perc_tab_deg()
    
    nums <- nums_tab_deg[c(1, 3, 2, 4)]
    percs <- perc_tab_deg[c(1, 3, 2, 4)]
    
    percs_lbls <- paste0(percs, "%")
    
    lab1 <- paste0("Mito_Down: ", nums[1])
    lab2 <- paste0("Mito_Up: ", nums[2])
    lab3 <- paste0("No-mito_Down: ", nums[3])
    lab4 <- paste0("No-mito_Up: ", nums[4])
    lbls <- c(lab1, lab2, lab3, lab4)
    
    main_pie <- paste("YOUR DATA\nTotal proteins:", sum(nums))
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(percs, radius = 0.8, labels = percs_lbls, col = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20))
  })
  
  # Drawing the pie plot with the proportions of mito and no_mito proteins AND the proportions of up and down regulated genes, where the total are differentially expressed ones. (download)
  pieDeg2 <- function(){
    nums_tab_deg <- nums_tab_deg()
    perc_tab_deg <- perc_tab_deg()
    
    nums <- nums_tab_deg[c(1, 3, 2, 4)]
    percs <- perc_tab_deg[c(1, 3, 2, 4)]
    
    percs_lbls <- paste0(percs, "%")
    
    lab1 <- paste0("Mito_Down: ", nums[1])
    lab2 <- paste0("Mito_Up: ", nums[2])
    lab3 <- paste0("No-mito_Down: ", nums[3])
    lab4 <- paste0("No-mito_Up: ", nums[4])
    lbls <- c(lab1, lab2, lab3, lab4)
    
    main_pie <- paste("YOUR DATA\nTotal proteins:", sum(nums))
    
    par(mai = c(0.62, 0.42, 0.42, 0.02))
    pie(percs, radius = 0.8, labels = percs_lbls, col = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20))
    title(main_pie, adj = 0.5, line = -1.2)
    legend("bottomleft", legend = lbls, fill = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20))
  }
  
  # Drawing a barplot comparing the percentages of the mito and no_mito proportions between the differentially expressed proteins and the whole proteome.
  bar_plot <- reactive({
    perc_new_mito <- perc_new_mito()
    perc_total_mito <- perc_total_mito()
    
    percents <- cbind(perc_new_mito, perc_total_mito)
    
    par(mai = c(0.82, 0.82, 0.82, 0.82))
    barplot(percents,
            names.arg = c("Your data", "The whole proteome"),
            main = "Comparing your data\n with the whole proteome",
            col = c("red", "cyan"),
            beside = T,
            ylim = c(0, 120))
    abline(h = 100, lty = 3)
    legend(3, 120, legend = c("Predicted Mitochondrial Location", "Predicted Ex-Mitochondrial Location"), fill = c("red", "cyan"), cex = 0.7)
  })
  
  # Drawing a barplot comparing the percentages of the mito and no_mito proportions between the differentially expressed proteins and the whole proteome. (download)
  bar_plot2 <- function(){
    perc_new_mito <- perc_new_mito()
    perc_total_mito <- perc_total_mito()
    
    percents <- cbind(perc_new_mito, perc_total_mito)
    
    par(mai = c(0.82, 0.82, 0.82, 0.82))
    barplot(percents,
            names.arg = c("Your data", "The whole proteome"),
            main = "Comparing your data\n with the whole proteome",
            col = c("red", "cyan"),
            beside = T,
            ylim = c(0, 120))
    abline(h = 100, lty = 3)
    legend(3, 120, legend = c("Mito", "No mito"), fill = c("red", "cyan"), cex = 0.7)
  }
  
  # Preparing the data to draw a GO plot.
  dataGO <- reactive({
    DatasetGO <- DatasetGO()
    
    dataGO <- matrix(t(DatasetGO[, 4:7]), nrow = 4)
    dataGO <- apply(dataGO, 2, as.numeric)
    colnames(dataGO) <- DatasetGO[, 1]
    rownames(dataGO) <- colnames(DatasetGO[, 4:7])
    dataGO
  })
  
  # Calculating the percentages of the GO plot.
  percsGO <- reactive({
    DatasetGO <- DatasetGO()
    dataGO <- dataGO()
    
    percsGO <- data.frame(numeric(nrow(DatasetGO)),
                          numeric(nrow(DatasetGO)),
                          numeric(nrow(DatasetGO)),
                          numeric(nrow(DatasetGO)),
                          stringsAsFactors = FALSE)
    colnames(percsGO) <- rownames(dataGO)
    for(i in 1:ncol(dataGO))
    {
      percsGO[i, ] <- round(dataGO[, i]/as.numeric(DatasetGO$TOTAL[i]) * 100, 2)
    }
    percsGO <- t(percsGO)
    colnames(percsGO) <- DatasetGO[, 1]
    percsGO
  })
  
  # Drawing a barplot with the percentages of the enriched/depleted GO terms.
  barplotGO <- reactive({
    percsGO <- percsGO()
    
    lbls <- c("Mito_Down", "Mito_Up", "No-mito_Down", "No-mito_Up")
    
    labGO <- NULL
    for(i in 1:ncol(percsGO))
    {
      name <- substr(colnames(percsGO)[i], 1, nchar(colnames(percsGO)[i]) - 13)
      code <- substr(colnames(percsGO)[i], nchar(colnames(percsGO)[i]) - 10, nchar(colnames(percsGO)[i]) - 1)
      if(nchar(name) < 50) labGO[i] <- paste(name, "\n", code, "")
      else if(nchar(name) < 100) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, nchar(name)), "\n", code, "")
      else if(nchar(name) < 150) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, nchar(name)), "\n", code, "")
      else if(nchar(name) < 200) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, 149), "\n", substr(name, 150, nchar(name)), "\n", code, "")
      else labGO[i] <- labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, 149), "\n", substr(name, 150, 199), "\n", substr(name, 200, nchar(name)), "\n", code, "")
    }
    
    tit <- paste("Top 15 GO terms in your dataset, and gene composition (in %)", "\n", "in terms of predicted cellular localization")
    
    par(mai = c(0.62, 3, 0.6, 0.08))
    yy <- barplot(percsGO[, 15:1], horiz = TRUE, names.arg = rep("", 15), col = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20), xlim = c(0, 115))
    title(tit, adj = 0.5, line = 0)
    text(x = 0, y = yy[1], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[15], cex = 0.8, pos = 2)
    text(x = 0, y = yy[2], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[14], cex = 0.8, pos = 2)
    text(x = 0, y = yy[3], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[13], cex = 0.8, pos = 2)
    text(x = 0, y = yy[4], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[12], cex = 0.8, pos = 2)
    text(x = 0, y = yy[5], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[11], cex = 0.8, pos = 2)
    text(x = 0, y = yy[6], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[10], cex = 0.8, pos = 2)
    text(x = 0, y = yy[7], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[9], cex = 0.8, pos = 2)
    text(x = 0, y = yy[8], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[8], cex = 0.8, pos = 2)
    text(x = 0, y = yy[9], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[7], cex = 0.8, pos = 2)
    text(x = 0, y = yy[10], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[6], cex = 0.8, pos = 2)
    text(x = 0, y = yy[11], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[5], cex = 0.8, pos = 2)
    text(x = 0, y = yy[12], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[4], cex = 0.8, pos = 2)
    text(x = 0, y = yy[13], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[3], cex = 0.8, pos = 2)
    text(x = 0, y = yy[14], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[2], cex = 0.8, pos = 2)
    text(x = 0, y = yy[15], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[1], cex = 0.8, pos = 2)
    text(x = coord_percs(percsGO[, 15]), y = yy[1], labels = lab_percs(percsGO[, 15]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 14]), y = yy[2], labels = lab_percs(percsGO[, 14]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 13]), y = yy[3], labels = lab_percs(percsGO[, 13]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 12]), y = yy[4], labels = lab_percs(percsGO[, 12]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 11]), y = yy[5], labels = lab_percs(percsGO[, 11]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 10]), y = yy[6], labels = lab_percs(percsGO[, 10]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 9]), y = yy[7], labels = lab_percs(percsGO[, 9]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 8]), y = yy[8], labels = lab_percs(percsGO[, 8]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 7]), y = yy[9], labels = lab_percs(percsGO[, 7]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 6]), y = yy[10], labels = lab_percs(percsGO[, 6]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 5]), y = yy[11], labels = lab_percs(percsGO[, 5]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 4]), y = yy[12], labels = lab_percs(percsGO[, 4]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 3]), y = yy[13], labels = lab_percs(percsGO[, 3]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 2]), y = yy[14], labels = lab_percs(percsGO[, 2]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 1]), y = yy[15], labels = lab_percs(percsGO[, 1]), pos = 2, cex = 0.7)
    legend(102, 6.8, legend = lbls, fill = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20), cex = 0.8)
  })
  
  # Drawing a barplot with the percentages of the enriched/depleted GO terms. (download)
  barplotGO2 <- function()
  {
    percsGO <- percsGO()
    
    lbls <- c("Mito_Down", "Mito_Up", "No-mito_Down", "No-mito_Up")
    
    labGO <- NULL
    for(i in 1:ncol(percsGO))
    {
      name <- substr(colnames(percsGO)[i], 1, nchar(colnames(percsGO)[i]) - 13)
      code <- substr(colnames(percsGO)[i], nchar(colnames(percsGO)[i]) - 10, nchar(colnames(percsGO)[i]) - 1)
      if(nchar(name) < 50) labGO[i] <- paste(name, "\n", code, "")
      else if(nchar(name) < 100) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, nchar(name)), "\n", code, "")
      else if(nchar(name) < 150) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, nchar(name)), "\n", code, "")
      else if(nchar(name) < 200) labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, 149), "\n", substr(name, 150, nchar(name)), "\n", code, "")
      else labGO[i] <- labGO[i] <- paste(substr(name, 1, 49), "\n", substr(name, 50, 99), "\n", substr(name, 100, 149), "\n", substr(name, 150, 199), "\n", substr(name, 200, nchar(name)), "\n", code, "")
    }
    
    tit <- paste("Top 15 GO terms in your dataset, and gene composition (in %)", "\n", "in terms of predicted cellular localization")
    
    par(mai = c(0.62, 3, 0.6, 0.08))
    yy <- barplot(percsGO[, 15:1], horiz = TRUE, names.arg = rep("", 15), col = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20), xlim = c(0, 115))
    title(tit, adj = 0.5, line = 0)
    text(x = 0, y = yy[1], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[15], cex = 0.8, pos = 2)
    text(x = 0, y = yy[2], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[14], cex = 0.8, pos = 2)
    text(x = 0, y = yy[3], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[13], cex = 0.8, pos = 2)
    text(x = 0, y = yy[4], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[12], cex = 0.8, pos = 2)
    text(x = 0, y = yy[5], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[11], cex = 0.8, pos = 2)
    text(x = 0, y = yy[6], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[10], cex = 0.8, pos = 2)
    text(x = 0, y = yy[7], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[9], cex = 0.8, pos = 2)
    text(x = 0, y = yy[8], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[8], cex = 0.8, pos = 2)
    text(x = 0, y = yy[9], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[7], cex = 0.8, pos = 2)
    text(x = 0, y = yy[10], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[6], cex = 0.8, pos = 2)
    text(x = 0, y = yy[11], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[5], cex = 0.8, pos = 2)
    text(x = 0, y = yy[12], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[4], cex = 0.8, pos = 2)
    text(x = 0, y = yy[13], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[3], cex = 0.8, pos = 2)
    text(x = 0, y = yy[14], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[2], cex = 0.8, pos = 2)
    text(x = 0, y = yy[15], srt = 0, adj = 1.1, xpd = TRUE, labels = labGO[1], cex = 0.8, pos = 2)
    text(x = coord_percs(percsGO[, 15]), y = yy[1], labels = lab_percs(percsGO[, 15]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 14]), y = yy[2], labels = lab_percs(percsGO[, 14]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 13]), y = yy[3], labels = lab_percs(percsGO[, 13]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 12]), y = yy[4], labels = lab_percs(percsGO[, 12]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 11]), y = yy[5], labels = lab_percs(percsGO[, 11]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 10]), y = yy[6], labels = lab_percs(percsGO[, 10]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 9]), y = yy[7], labels = lab_percs(percsGO[, 9]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 8]), y = yy[8], labels = lab_percs(percsGO[, 8]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 7]), y = yy[9], labels = lab_percs(percsGO[, 7]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 6]), y = yy[10], labels = lab_percs(percsGO[, 6]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 5]), y = yy[11], labels = lab_percs(percsGO[, 5]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 4]), y = yy[12], labels = lab_percs(percsGO[, 4]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 3]), y = yy[13], labels = lab_percs(percsGO[, 3]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 2]), y = yy[14], labels = lab_percs(percsGO[, 2]), pos = 2, cex = 0.7)
    text(x = coord_percs(percsGO[, 1]), y = yy[15], labels = lab_percs(percsGO[, 1]), pos = 2, cex = 0.7)
    legend(102, 6.8, legend = lbls, fill = c("red", "red", "cyan", "cyan"), density = c(NA, 20, NA, 20), cex = 0.8)
  }
  
###### OUTPUTS ######
  
  # Taking the column names of the input data to select one of them as the Log Fold Change column.
  observe({
    req(input$microarray)
    dsnames <- names(csv())
    cb_options <- list()
    cb_options[dsnames] <- dsnames
    updateSelectInput(session, "cols_logfc",
                      label = "Select column with Log Fold Change data:",
                      choices = cb_options,
                      selected = "")
  })
  
  # Taking the column names of the input data to select one of them as the p-value column.
  observe({
    req(input$microarray)
    dsnames <- names(csv())
    cb_options <- list()
    cb_options[dsnames] <- dsnames
    updateSelectInput(session, "cols_pval",
                      label = "Select column with p-value data:",
                      choices = cb_options,
                      selected = "")
  })
  
  # Taking the column names of the input data to select one of them as the Entrez Gene ID column.
  observe({
    req(input$microarray)
    dsnames <- names(csv())
    cb_options <- list()
    cb_options[dsnames] <- dsnames
    updateSelectInput(session, "cols_entrez",
                      label = "Select column with Entrez Gene ID data:",
                      choices = cb_options,
                      selected = "")
  })

#### Dimensions ####
  
  # Initial data
  output$DimInitial <- renderText({
    if(!is.null(input$microarray))
      paste("There are", dim(csv())[1], "rows in the initial table.")
  })
  
  # Proteins up & down regulated, simultaneously, data
  output$DimUpDown <- renderText({
    if(!is.null(input$microarray))
      paste(dim(up_down_dataframe())[1], "proteins in your dataset are found to be negatively and positively expressed, simultaneously. They have been removed the analysis, but can be downloaded")
  })
  
  # Differentually expressed genes data
  output$DimRem_Use <- renderText({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      paste("We have found", dim(rem_use())[1], "differentially expressed genes.")
  })
  
  # Differentually expressed proteins data
  output$DimFinalData <- renderText({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      paste("We have found", dim(mega_dataframe())[1], "differentially expressed proteins.")
  })
  
  # Gene Ontology data
  output$DimDFGO <- renderText({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      paste("There are", dim(DatasetGO())[1], "Gene Ontology terms affected by the differentially expressed genes.")
  })

#### Dataframes

  # Initial data
  output$Input <- renderDataTable({
      if(!is.null(input$microarray))
        csv()
    })
  
  # Differentually expressed genes data
  output$Rem_Use <- renderDataTable({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      rem_use()
  })
  
  # Differentually expressed proteins data
  output$FinalData <- renderDataTable({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      mega_dataframe()
  })
  
  # Gene Ontology data
  output$DFGO <- renderDataTable({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      DatasetGO()
  })

#### Conclusions ####
  
  # p-value of the mito proportion test.
  output$PValTest <- renderText({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      conclusions()[[1]]
  })
  
  # Conclusion phrase of the mito proportion test.
  output$Conclusion <- renderText({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      conclusions()[[2]]
  })

#### Graphs ####
  
  # Volcano plot
  output$Volcano <- renderPlot({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      volc()
  })
  
  # Pie plot with the proportions of mito and no_mito proteins, where the total are the differentially expressed ones.
  output$PieNew <- renderPlot({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      pieNew()
    })
  
  # Pie plot with the proportions of mito and no_mito proteins, where the total are the whole proteome ones.
  output$PieTotal <- renderPlot({
    if(!is.null(input$microarray))
       pieTotal()
    })
  
  # Pie plot with the proportions of mito and no_mito proteins AND the proportions of up and down regulated genes, where the total are differentially expressed ones.
  output$PieDeg <- renderPlot({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      pieDeg()
  })
  
  # Barplot comparing the percentages of the mito and no_mito proportions between the differentially expressed proteins and the whole proteome.
  output$Barplot <- renderPlot({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      bar_plot()
    })
  
  # Barplot with the percentages of the enriched/depleted GO terms.
  output$BarplotGO <- renderPlot({
    if((!is.null(input$microarray)) && (!is.null(input$cols_logfc)) && (!is.null(input$cols_pval)) && (!is.null(input$cols_entrez)))
      barplotGO()
  })

#### Downloads ####
  
  # Volcano plot (download)
  output$VolcanoDown <- downloadHandler(
    filename = "VolcanoPlot.png",
    content = function(file) {
      png(file)
      volc2()
      dev.off()
    },
    contentType = 'png')

  # Pie plot with the proportions of mito and no_mito proteins, where the total are the differentially expressed ones. (download)
  output$PieNewDown <- downloadHandler(
    filename = "Your_proteins_PIE.png",
    content = function(file) {
      png(file)
      pieNew2()
      dev.off()
    },
    contentType = 'png')
  
  # Pie plot with the proportions of mito and no_mito proteins, where the total are the whole proteome ones. (download)
  output$PieTotalDown <- downloadHandler(
    filename = "The_whole_proteome_PIE.png",
    content = function(file) {
      png(file)
      pieTotal2()
      dev.off()
    },
    contentType = 'png')
  
  # Pie plot with the proportions of mito and no_mito proteins AND the proportions of up and down regulated genes, where the total are differentially expressed ones. (download)
  output$PieDegDown <- downloadHandler(
    filename = "Targeting_DEG_PIE.png",
    content = function(file) {
      png(file)
      pieDeg2()
      dev.off()
    },
    contentType = 'png')
  
  # Barplot comparing the percentages of the mito and no_mito proportions between the differentially expressed proteins and the whole proteome. (download)
  output$BarplotDown <- downloadHandler(
    filename = "Comparison_Barplot.png",
    content = function(file) {
      png(file)
      bar_plot2()
      dev.off()
    },
    contentType = 'png')
  
  # Barplot with the percentages of the enriched/depleted GO terms. (download)
  output$downBarplotGO <- downloadHandler(
    filename = "GO_Barplot.png",
    content = function(file) {
      png(file)
      barplotGO2()
      dev.off()
    },
    contentType = 'png')
  
  # Differentually expressed genes data (download)
  output$downRem_Use <- downloadHandler(
    filename = "Differentially_expressed_genes.csv",
    content = function(file) {
      write.csv(rem_use(), file, row.names = FALSE)
    }
  )
  
  # Differentually expressed proteins data (download)
  output$downFinalData <- downloadHandler(
    filename = "Final_data.csv",
    content = function(file) {
      write.csv(mega_dataframe(), file, row.names = FALSE)
    }
  )
  
  # Proteins up & down regulated, simultaneously, data (download)
  output$downUpDownData <- downloadHandler(
    filename = "Genes_up_and_down_reg.csv",
    content = function(file) {
      write.csv(up_down_dataframe(), file, row.names = FALSE)
    }
  )
  
  # Gene Ontology data (download)
  output$downDFGO <- downloadHandler(
    filename = "GO_table.csv",
    content = function(file) {
      write.csv(DatasetGO(), file, row.names = FALSE)
    }
  )
}