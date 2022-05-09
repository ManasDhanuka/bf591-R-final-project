options(shiny.maxRequestSize=30*1024^2)
library(shiny)
library(ggplot2)
library(colourpicker)
library(dplyr)
library(shinythemes)
library(DT)
library(tidyverse)
library(data.table)
library(matrixStats)
library(gplots)
library(matrixTests)
library(tidyr)
library(fgsea)
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library(Rcpp)
library(ComplexHeatmap)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),
    
    # Application title
    titlePanel("BF591 Final Project"),
    
    ###########################USER INTERFACE###################################
    
    tabsetPanel(
        ###############SUMMARY SECTION UI################
        tabPanel(title = "Summary", fluid = TRUE, br(),
                 sidebarLayout(
                     sidebarPanel( width = 3,fileInput(inputId = "summary_datafile", "Upload the summary file", accept = ".csv", 
                                                       buttonLabel = "Browse", placeholder = "sample_info.csv")
                                
                                   
                                   
                     ),
                     
                     mainPanel(
                         tabsetPanel(
                             tabPanel(title = "Summary", br() ,tableOutput("summary_table")),
                             tabPanel(title = "Table", br(), DT::dataTableOutput("meta_data_table")),
                             tabPanel(title = "Plot", br(),
                                      
                                      sidebarLayout(
                                          sidebarPanel(
                                              width = 3,
                                              radioButtons(inputId = "meta_plot_para",label = "Select the value to be plotted as a histogram",
                                                           choices = c("age_of_onset","mrna.seq_reads","rin", "cag", "Duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade"),
                                                           selected = "age_of_onset")

                                              
                                          ),#plot selection panel
                                          mainPanel(plotOutput("summary_plot")) #summary plot output main panel
                                      )
                                      
                                      
                                      
                                      
                             )
                         )
                     ) #overall summary main panel
                     
                     
                     # closing sidebarPanel() 
                 ) ), # "Summary "tabPanel() closing
        
        ###############COUNTS MATRIX SECTION UI#################
        tabPanel(title = "Counts Matrix Analysis", br(), 
                 sidebarLayout(
                     
                     #input elements 
                     sidebarPanel( width = 3,
                                   
                                   #counts file input
                                   fileInput(inputId = "counts_datafile", "Upload the counts matrix file", accept = ".csv", 
                                             buttonLabel = "Browse", placeholder = "norm_counts.csv"),
                                   
                                   #slider inputs
                                   sliderInput(inputId = "var_slider", label = "Adjust the percentile variance to be included", min = 0, max = 100, value = 50),
                                   sliderInput(inputId = "non_zero_slider", label = "Adjust the level non-zero genes present" , min = 0, max = 50, value = 20)
                                   
                                   
                                   
                                   
                                   
                     ),
                     
                     mainPanel(
                         tabsetPanel(
                             tabPanel(title = "Summary", br(),
                                      mainPanel(tableOutput("cts_summary_table"))),
                             tabPanel(title = "Diagnostic Plots", br(),
                                      tabsetPanel(
                                          tabPanel(title = "Median Count vs Variance", br(),
                                                   mainPanel(plotOutput("dp_cts_1"), width = "100%")),
                                          tabPanel(title = "Median Count vs Number of Zeros", br(), 
                                                   mainPanel(plotOutput("dp_cts_2"), width = "100%"))
                                      )
                             ),
                             tabPanel(title = "Heatmap", br(),
                                      
                                      sidebarLayout(
                                          
                                          sidebarPanel(
                                              
                                              radioButtons(inputId = "log_transform", "Apply Log Transformation to Counts",
                                                           choices = c("Yes", "No"), selected = "Yes")
                                              
                                              
                                          ),
                                          mainPanel(plotOutput("cts_heatmap"))
                                      )
                             ),
                             
                             tabPanel(title = "PCA", br(),
                                      
                                      sidebarLayout(
                                          sidebarPanel( width = 3,
                                                        
                                                        radioButtons(inputId = "first_PC", "Select the PC on X-axis",
                                                                     choices = c("PC1","PC2","PC3","PC4","PC5"),
                                                                     selected = "PC1"),
                                                        
                                                        radioButtons(inputId = "second_PC", "Select the PC on Y-axis",
                                                                     choices = c("PC1","PC2","PC3","PC4","PC5"),
                                                                     selected = "PC2")
                                                        
                                                        
                                                        
                                          ),
                                          
                                          mainPanel(plotOutput("pca_plot"))
                                      )
                                      
                             )
                         )
                     )
                     
                     
                     # closing sidebarPanel() 
                 ) ), # "Summary "tabPanel() closing),
        
        ############### DIFFERENTIAL EXPRESSION UI #############
        tabPanel(title = "Differential Expression", br(),
                 
                 sidebarLayout(
                     
                     sidebarPanel( width = 3,
                                   
                                   #taking file input 
                                   fileInput(inputId = "de_datafile", "Upload a Differential Gene Expession file", 
                                             accept = ".csv", 
                                             buttonLabel = "Browse", placeholder = "deseq_res.csv"),
                                   
                                   h5("A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis."),
                                   
                                   
                                   #values for the x and y axis 
                                   #x value
                                   radioButtons(inputId = "x_axis",label = "Select the value to be plotted on the x axis",
                                                choices = c("baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj"),
                                                selected = "log2FoldChange"),
                                   
                                   #y value
                                   radioButtons(inputId = "y_axis",label = "Select the value to be plotted on the y axis",
                                                choices = c("baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj"),
                                                selected = "padj"),
                                   #Color pickers
                                   #1
                                   colourInput(
                                       inputId = "base_point_color",
                                       labe = "Pick a base point color:",
                                       value = "#226764",
                                       showColour = c("both","text","background"),
                                       palette = c("square"),
                                       allowedCols = NULL,
                                       allowTransparent = FALSE,
                                       returnName = FALSE,
                                       closeOnClick = FALSE),
                                   
                                   #2
                                   colourInput(
                                       inputId = "high_point_color",
                                       labe = "Pick a high point color:",
                                       value = "#A8383B",
                                       showColour = c("both","text","background"),
                                       palette = c("square"),
                                       allowedCols = NULL,
                                       allowTransparent = FALSE,
                                       returnName = FALSE,
                                       closeOnClick = FALSE),
                                   
                                   #pvalue selecter
                                   sliderInput(inputId= "pvalue_slider", label= "Select the magnitude of the p adjusted coloring:", 
                                               min = -35, max= 0, value = -4)
                                   #Submit button
                                   
                     ),
                     
                     mainPanel(
                         
                         tabsetPanel(
                             tabPanel("Plot", br(), plotOutput("volcano"), width = "100%"), 
                             tabPanel("Table", br(), DT::dataTableOutput("de_table"))
                             
                         )#closing the tabsetPanel()
                         
                     )#closing the mainPanel()
                     
                     
                 )),
        
        ###################### GSEA UI #####################
        
        tabPanel(title = "GSEA", br(),
                 
                 sidebarLayout(
                     sidebarPanel( width = 3,
                                   
                                   fileInput(inputId = "gsea_de_datafile", "Upload a Differential Gene Expession file", 
                                             accept = ".csv", 
                                             buttonLabel = "Browse", placeholder = "deseq_res.csv")
                                   
                     ),
                     mainPanel(
                         tabsetPanel(
                             
                             tabPanel("NES Bar Plot", br(), 
                                      sidebarLayout(
                                          sidebarPanel( width = 4,
                                                        sliderInput(inputId = "gsea_nop_slider", label = "Select the number of Top Pathways to be plotted",
                                                                    min = 1, max = 20, value = 10)
                                                        
                                                        
                                                        
                                          ),
                                          mainPanel(plotOutput("pathways_barplot")) #output horizontal bar plot of gsea
                                      ) 
                                      
                             ),
                             tabPanel("Table", br(), 
                                      
                                      sidebarLayout(
                                          sidebarPanel(width = 4,
                                                       
                                                       sliderInput(inputId = "gsea_padj_slider", label = "Select the padj threshold",
                                                                   min = -22, max = 0, value = -10),
                                                       
                                                       radioButtons(inputId = "pathway_type", label = "Select the NES pathway type",
                                                                    choices = c("All Pathways","Positive NES Pathways Only","Negative NES Pathways Only"),
                                                                    selected = "All Pathways"),
                                                       
                                                       
                                                       
                                                       downloadButton(outputId ="download_gsea_table",label = "Download", width = "10%")
                                          ), 
                                          mainPanel(DT::dataTableOutput("gsea_filtered_table"))
                                      )
                                      
                                      
                             ),
                             tabPanel("NES Scatter Plot", br(), 
                                      
                                      sidebarLayout(
                                          sidebarPanel(
                                              sliderInput(inputId = "gsea_sp_slider", label = "Select the padj threshold",
                                                          min = -22, max = 0, value = -10)
                                              
                                             
                                              
                                          ),
                                          
                                          mainPanel(plotOutput("gsea_scatter_plot"))
                                      )
                             )
                         )
                     )
                     
                 ))
    ) #tabsetPanel closing
    
)#closing fluidpage()

server <- function(input, output) {
    
    ###############################FUNCTIONS####################################
    
    ######SUMMARY FUNCTIONS#######
    
    #Loading meta data
    load_meta_data <- reactive({
        
        req(input$summary_datafile)
        inFile <- input$summary_datafile
        if (is.null(inFile)){
            return(NULL)  }
        else{
            meta_data_df <- read.csv(inFile$datapath, header=TRUE, sep=",")
        }
        return(meta_data_df) 
        
    })
    
    #function to produce summary table 
    summary_stats <- function(meta_df){
        
        meta_summ <- dplyr::select(meta_df, c("age_of_death", "Diagnosis", "Organism", "tissue", "age_of_onset", "cag", "Duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade", "mrna.seq_reads","rin"))
        
        Column_Name <- c("Diagnosis","tissue","age_of_onset","mrna.seq_reads","rin", "cag", "Duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade")
        
        Data_type <- c(typeof(meta_summ$Diagnosis[1]), typeof(meta_summ$tissue[1]), typeof(meta_summ$age_of_onset[1]), typeof(meta_summ$mrna.seq_reads[1]), typeof(meta_summ$rin[1]),typeof(meta_summ$cag[1]),typeof(meta_summ$Duration[1]),typeof(meta_summ$h.v_cortical_score[1]),typeof(meta_summ$h.v_striatal_score[1]),typeof(meta_summ$vonsattel_grade[1]))
        
        Mean_or_Distinct_values <- c(paste(unique(meta_summ$Diagnosis)[1],"or",unique(meta_summ$Diagnosis)[2]),c(unique(meta_summ$tissue)),mean(meta_summ$age_of_onset, na.rm = TRUE),mean(meta_summ$mrna.seq_reads, na.rm = TRUE),mean(meta_summ$rin, na.rm = TRUE),mean(meta_summ$cag, na.rm = TRUE),mean(meta_summ$Duration, na.rm = TRUE),mean(meta_summ$h.v_cortical_score, na.rm = TRUE),mean(meta_summ$h.v_striatal_score, na.rm = TRUE),paste(unique(meta_summ$vonsattel_grade)[2],"or",unique(meta_summ$vonsattel_grade)[3]))
        
        summary_df <- data.frame(Column_Name,Data_type,Mean_or_Distinct_values)
        
        colnames(summary_df) <- c("Column Name", "Data Type", "Mean or Distinct Values")
        
        return(summary_df)
    }
    
    #function to produce meta data plots
    meta_data_plot<- function(meta_df,x_para){
        
        meta_plot <- ggplot(meta_df, aes(x=!!sym(x_para)))+
            geom_histogram(color="#0D4D4B", fill="#679B99") + theme_light(base_size = 16)
        
        return(meta_plot)
    }
    
    ########COUNTS FUNCTIONS#########
    
    #loading counts data
    load_counts_data <- reactive({
        
        req(input$counts_datafile)
        inFile <- input$counts_datafile
        if (is.null(inFile)){
            return(NULL)  }
        else{
            cts_data_df <- read.csv(inFile$datapath, header=TRUE, sep="\t")
            rownames(cts_data_df) <- cts_data_df[,1]
            cts_data_df <- cts_data_df[,-1]
        }
        return(cts_data_df) 
        
    })
    
    #function to filter data depending upon the parametrs selected 
    data_filter <- function(cts_data, cutoff, zeros){
        
        cutoff <- cutoff/100
        variances <- apply(cts_data, 1, var)
        cutoff <- quantile(variances, probs = cutoff)
        store <- apply(cts_data, 1, function(x){var(x)>cutoff&sum(x==0)<zeros})
        fil_cts_data <- cts_data[store,]
        
        return(fil_cts_data)
    }
    
    #function to create summarized result after filtering
    cts_summ_table <- function(cts_data, cutoff, zeros){
        
        #filtering cts_data
        filtered_cts_data <- data_filter(cts_data = cts_data, cutoff, zeros) 
        
        #number of samples
        num_samples <- ncol(cts_data)
        
        #total number of genes
        num_genes <- nrow(cts_data)
        
        # number and % of genes passing current filter
        num_pass_fil <- nrow(filtered_cts_data)
        per_pass_fil <- (num_pass_fil/num_genes)*100
        
        # number and % of genes not passing current filter
        not_pass_fil <- num_genes - num_pass_fil
        per_not_pass_fil <- (not_pass_fil/num_genes)*100
        
        Conditon <- c("number of samples", "total number of genes", 
                      "number of genes passing current filter", "percentage of genes passsing current filter",
                      "number of genes not passing current filter", "percentage of genes not passing the current filter")
        
        Count <- c(num_samples,num_genes, num_pass_fil, capture.output(cat(per_pass_fil,'%')), not_pass_fil, capture.output(cat(per_not_pass_fil,'%')) )
        
        filter_summ_df <- data.frame(Conditon,Count)
        
        return(filter_summ_df)
        
    }
    
    #function to return scatter plots of median counts vs variance and median counts vs number of zeros
    diagnostic_plot <- function(data, cuttoff, zeros_cuttoff){
        
        #calculating diagnostic parameters
        cuttoff <- cuttoff/100
        variances <- apply(data, 1, var)
        cuttoff <- quantile(variances, probs = cuttoff)
        zeros <- apply(data, 1, function(x){sum(x==0)})
        medians <- apply(data, 1, median)
        
        plotData <- tibble(variance = variances,
                           zero = zeros,
                           median = medians)
        
        # plotting data
        var_plot <- plotData %>%
            ggplot(aes(x = log(median), y = log(variance), colour = variance<cuttoff)) +
            geom_point(stat = "identity") +
            scale_colour_manual(values = c("#003432", "#679B99")) +
            ylab("Variance") + xlab("log(Median)") +
            ggtitle("Log(Median Normalized Counts) vs Log(Variance Normalized Counts)") +
            theme(legend.position="bottom") + theme_light(base_size = 16)
        
        zero_plot <- plotData %>%
            ggplot(aes(x = log(median), y = zero, colour = zero>zeros_cuttoff)) +
            geom_point(stat = "identity") +
            scale_colour_manual(values = c("#003432", "#679B99")) +
            ylab("Number of Zeros") + xlab("log(Median)") +
            ggtitle("Log(Median Normalized Counts) vs Number of Zero Values") +
            theme(legend.position="bottom") + theme_light(base_size = 16)
            
        
        return(list(var_plot, zero_plot))
    }
    
    #function to plot heatmap
    # plot_heatmap <- function(cts_data, cutoff, zeros){
    # 
    #     filtered_data <- data_filter(cts_data, cutoff, zeros)
    # 
    #     #code block to color clusters
    #     s <- colnames(cts_data)
    #     type = append(rep(c("Neurologically Normal"), 49), rep(c("Huntingtons Disease"), 20))
    #     condition_colors<- data.frame(s,type)
    #     condition_colors <- transmute(condition_colors, color = if_else(type == "Neurologically Normal", "blue", "red"))
    # 
    #     hm_data <- data.matrix(filtered_data)
    # 
    #     #plotting heatmap
    #     hm <- heatmap.2(hm_data,ColSideColors=condition_colors$color, xlab="Samples", ylab="Genes",
    #                     main="Gene Expression Across Samples",trace="none", density.info = "none",
    #                     key.xlab="Expression Level", scale="row", cexRow = 0.5, cexCol = 0.2, key = TRUE)
    #          legend(x=0.9, y=1, xpd=TRUE, inset=c(-0.15,0),
    #            legend=c("Huntington's Disease", 'Neurologically Normal'), title='Sample Type', fill=c('red', 'blue'))
    # 
    #     return(hm)
    # }
    # 
    
    plot_heatmap <- function(cts_data, cutoff , zeros, log_transform = "Yes"){
        
        
        filtered_df <- data_filter(cts_data, cutoff, zeros)
        
        hp_data <- as.matrix(filtered_df)
        
        if(log_transform == "No"){
            
            hplot <- Heatmap(hp_data, name = 'Counts', cluster_rows = FALSE,
                             column_gap = unit(0, "mm"), show_row_names = FALSE, border = TRUE)
            
        }else{
            
            # log_plot_data <- log10(hp_data+1)
            
            hplot <- Heatmap(log(hp_data+1), name = 'log10(Counts+1)', cluster_rows = FALSE,
                             column_gap = unit(0, "mm"), show_row_names = FALSE, border = TRUE)
            
        }
        return(hplot)
    }
    
    #function to calculate and produce PCA plot
    
    plot_PCA <- function(cts_data,cuttoff, zeros, x_pca, y_pca){
        
        #filtering data
        filtered_data <- data_filter(cts_data, cuttoff, zeros)
        
        #calculate principal components
        pca_results <- prcomp(filtered_data)
        
        pca_df <- as.data.frame(pca_results$x)
        
        #calculate percent variance
        pca_vars <- pca_results$sdev^2 / sum(pca_results$sdev^2)
        
        #creating labels for axes
        x_var <- paste(c(round((pca_vars[match(x_pca,colnames(pca_df))])*100,2),"%"),collapse = "")
        y_var <- paste(c(round((pca_vars[match(y_pca,colnames(pca_df))])*100,2),"%"),collapse = "")
        
        x_label <- paste(c(x_pca,"(", x_var, ")"), collapse = " ")
        y_label <- paste(c(y_pca,"(", y_var, ")"), collapse = " ")
        
        
        pca_plot <- ggplot(pca_df, aes(x = !!sym(x_pca), y= !!sym(y_pca))) + geom_point(size = 3, color = "#41817F")+
            ggtitle(paste(c("Plot of", x_pca, "vs", y_pca), collapse = " "))+
            xlab(x_label) + ylab(y_label) + theme_light(base_size = 16)
        
        return(pca_plot)
    }
    
    #########DIFFERENTIAL EXPRESSION FUNCTIONS#########
    
    load_de_data <- reactive({
        
        req(input$de_datafile)
        inFile <- input$de_datafile
        if (is.null(inFile)){
            return(NULL)  }
        else{
            data_df <- read.csv(inFile$datapath, header=TRUE, sep="\t")
        }
        return(data_df) 
        
    })
    
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        
        plotting_data <- mutate(dataf, 
                                padj_thresh = if_else( padj > 1 * 10^slider, TRUE, FALSE))
        
        vol_plot <- ggplot(data = plotting_data, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), colour = padj_thresh)) + geom_point() + 
            scale_color_manual(name = paste("padj <", as.character(1 * 10^slider)),values = c(color1,color2) ) + theme_light(base_size=16)+
            theme(legend.position="bottom")
        
        return(vol_plot)
    }
    
    draw_table <- function(dataf, slider) {
        
        colnames(dataf)[1] <- "gene"
        
        data_df_filtered <- filter(dataf, dataf$padj < 1 * 10^slider )
        
        j1 <- c(7)
        data_df_filtered[j1] <- lapply(data_df_filtered[j1], formatC, format = "e", digits = 5)
        
        j2 <- c(3:6)
        data_df_filtered[j2] <- lapply(data_df_filtered[j2], formatC, format = "f", digits = 3)
        
        j3<- c(8:10)
        data_df_filtered[j3] <- lapply(data_df_filtered[j3], formatC, format = "f", digits = 3)
        
        return(data_df_filtered)
    }
    
    ########### GSEA FUNCTIONS #############
    
    load_gsea_de_data <- reactive({
        
        req(input$gsea_de_datafile)
        inFile <- input$gsea_de_datafile
        if (is.null(inFile)){
            return(NULL)  }
        else{
            data_df <- read.csv(inFile$datapath, header=TRUE, sep="\t")
        }
        return(data_df) 
        
    })
    
    #function to run fgsea() on differential expression data
    run_fgsea <- function(de_data){
        
        pathways.hallmark <- gmtPathways("Data/h.all.v7.5.1.symbols.gmt")
        
        ranks <- de_data %>%
            drop_na(symbol, log2FoldChange) %>%
            distinct(symbol, log2FoldChange, .keep_all=TRUE) %>%
            arrange(desc(log2FoldChange)) %>%
            dplyr::select(symbol, log2FoldChange) %>%
            deframe()
        
        fgsea_results <- as_tibble(fgsea(pathways.hallmark, ranks))
        
        return(fgsea_results)
        
    }
    
    #function to plot the selected number of pathways
    top_pathways <- function(de_data, num_paths){
        
        #perform fgsea on the expression data first
        fgsea_results <- run_fgsea(de_data)
        
        #Extracting top negative and positive NES scored pathways
        top_pos_NES <- fgsea_results %>% slice_max(NES, n=num_paths) %>% pull(pathway)
        top_neg_NES <- fgsea_results %>% slice_min(NES, n=num_paths) %>% pull(pathway)
        
        #Creating table containing pathways to plot (top positive and top negative)
        to_plot <- fgsea_results %>% 
            filter(pathway %in% c(top_pos_NES, top_neg_NES)) %>%
            mutate(pathway = factor(pathway)) %>%
            #Removing _ to plot the pathway name properly
            mutate(plot_name = str_replace_all(pathway, '_', ' '))
        
        plot <- to_plot %>% 
            #Reordering by NES
            mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
            ggplot() +
            geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
            scale_fill_manual(values = c('TRUE' = '#A8383B', 'FALSE' = '#226764')) + 
            theme_light(base_size = 12) +
            ggtitle('fgsea results for Hallmark MSigDB gene sets') +
            ylab('Normalized Enrichment Score (NES)') +
            xlab('') +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
            #to make graph horizontal
            coord_flip() 
        
        return(plot)
    }
    
    #function to filter table depending upon the padj value
    draw_gsea_table <- function(de_data, slider, NES_type = "All Pathways"){
        
        fgsea_results <- run_fgsea(de_data)
        
        if (NES_type == "Negative NES Pathways Only") {
            
            filtered_data <- filter(fgsea_results, fgsea_results$padj < 1 * 10^slider & fgsea_results$NES < 0)
            
        } else if (NES_type == "Positive NES Pathways Only" ) {
            
            filtered_data <- filter(fgsea_results, fgsea_results$padj < 1 * 10^slider & fgsea_results$NES > 0)
            
        } else {
            
            filtered_data <- filter(fgsea_results, fgsea_results$padj < 1 * 10^slider)
        }
        
        j1 <- c(2:4)
        filtered_data[j1] <- lapply(filtered_data[j1], formatC, format = "e", digits = 5)
        
        j2<- c(5:6)
        filtered_data[j2] <- lapply(filtered_data[j2], formatC, format = "f", digits = 3)
        
        
        
        
        return(filtered_data[1:6])
        
    }
    
    gsea_scatter_plot <- function(de_data, slider){
        
        fgsea_results <- run_fgsea(de_data)
        
        plotting_data <- mutate(fgsea_results, color = if_else(padj < 1*10^slider, 'grey', '#41817F' ))
        
        sp <- ggplot(plotting_data, aes(x = NES, y = -log10(pval), color = color)) + 
            scale_color_manual(name = paste("padj <", as.character(1*10^slider)),values = c("grey","#41817F"))+
            theme(legend.position="bottom") + geom_point(stat= "identity", show.legend = FALSE) + 
            ylab("NES Values") + xlab("-log10(pval)") +
            ggtitle("-Log10(pval) vs NES Values") + theme_light(base_size = 16)
        
        return(sp)
        
    }
    
    
    ###############################OUTPUTS######################################
    
    ####### SUMMARY OUTPUTS #######
    
    #call to summary table output
    output$summary_table <- renderTable({
        summary_dataf <- load_meta_data()
        summary_stats(summary_dataf)
    }) #
    
    #call to full meta_data table output
    output$meta_data_table <-DT::renderDataTable({
        load_meta_data()
    }) #
    
    #call to meta_data plot outout 
    output$summary_plot <- renderPlot({
        
        meta_df <- load_meta_data()
        
        meta_data_plot(meta_df,x_para = input$meta_plot_para )
        
    })
    
    ########## COUNTS OUTPUTS ################
    
    #filtered summary
    output$cts_summary_table <- renderTable({
        cts_summ_table(load_counts_data(), input$var_slider, input$non_zero_slider)
    })
    
    #median vs variance plot
    output$dp_cts_1 <- renderPlot({
        diagnostic_plot(load_counts_data(), input$var_slider, input$non_zero_slider)[1]
    })
    
    #median vs number of zeros
    output$dp_cts_2 <- renderPlot({
        diagnostic_plot(load_counts_data(), input$var_slider, input$non_zero_slider)[2]
    })
    
    #heatmap
    output$cts_heatmap <- renderPlot({
        plot_heatmap(load_counts_data(), input$var_slider, input$non_zero_slider, input$log_transform)
        # plot(x=10, y=10)
    }, height = 500, width = 800)
    
    #PCA plot
    output$pca_plot <- renderPlot({
        plot_PCA(load_counts_data(), input$var_slider, input$non_zero_slider, input$first_PC, input$second_PC)
    })
    
    ####### DIFFERENTIAL EXPRESSION OUTPUTS #########
    
    #differential expression volcano plot
    output$volcano <- renderPlot({
        volcano_plot(load_de_data(), x_name = input$x_axis, y_name = input$y_axis, slider = input$pvalue_slider,
                     color1 = input$base_point_color, color2= input$high_point_color)
        
    }, height = 600, width = 800)
    
    #differential expression table after filters
    output$de_table <- DT::renderDataTable({
        dataf <- load_de_data()
        draw_table(dataf, slider = input$pvalue_slider) 
    })
    
    ############### GSEA OUTPUTS ###############
    
    #NES horzontal barplot
    output$pathways_barplot <- renderPlot({
        top_pathways(load_gsea_de_data(), num_paths = input$gsea_nop_slider )
    })
    
    #filtered table
    output$gsea_filtered_table <- DT::renderDataTable({
        draw_gsea_table(load_gsea_de_data(), slider = input$gsea_padj_slider, NES_type = input$pathway_type )
    })
    
    #download handler
    output$download_gsea_table <- downloadHandler(
        filename = function(){"Results.csv"},
        content = function(filename){
            write.csv(draw_gsea_table(load_gsea_de_data(), slider = input$gsea_padj_slider, NES_type = input$pathway_type ), filename)
        }
    )
    
    output$gsea_scatter_plot <- renderPlot({
        gsea_scatter_plot(load_gsea_de_data(), slider = input$gsea_sp_slider)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
