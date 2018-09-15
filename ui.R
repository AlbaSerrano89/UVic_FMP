shinyUI(
  fluidPage(
    titlePanel("MicroGOTP"),
    conditionalPanel(
      condition = "input.sidebar == false",
      fluidRow(
        column(4,
               fileInput("microarray", "Input data (CSV file)", accept = c("text/csv", ".csv", "text/comma-separated-values,text/plain")),
               helpText("The file must contain, at least, these three columns:"),
               helpText("- A column with the Entrez Gene IDs"),
               helpText("- A column with the LogFC values of the microarray analysis"),
               helpText("- A column with the p-values of the microarray analysis"),
               br(),
               selectInput("specie", "Species",
                           choices = list("H. sapiens" = 9606,
                                          "M. musculus" = 10090,
                                          "A. thaliana" = 3702,
                                          "D. melanogaster" = 7227,
                                          "C. elegans" = 6239,
                                          "S. cerevisiae" = 559292
                                          ),
                           selected = 7227)
               ),
        column(4,
               h3("Expression"),
               selectInput("cols_logfc", "Select column with Log Fold Change data:", c("")),
               numericInput("logfc", label = "What Log Fold Change threshold would you like to use?", value = 1.5),
               selectInput("cols_pval", "Select column with p-value data:", c("")),
               numericInput("pval_init", label = "What significance threshold would you like to use?", value = 0.05),
               selectInput("cols_entrez", "Select column with Entrez Gene ID data:", c(""))
               ),
        column(4,
               h3("Redundant transcripts"),
               selectInput("numb_genes", "Redundant transcripts in data set?",
                           c("Remove all of them" = "rem", "Use all of them" = "use")),
               br(),
               h3("Confidence of the proportion tests"),
               numericInput("confidence", label = "For the localization of proteins", value = 95),
               br(),
               numericInput("confidenceGO", label = "For the GO analysis", value = 95)
               )
        ),
      fluidRow(
        h3("Which information do you want?"),
        conditionalPanel(
          condition = "input.specie == 9606",
          checkboxGroupInput("human_var", label = "",
                             choices = list("UCSC" = 3,
                                            "PDB" = 4,
                                            "Reactome" = 5,
                                            "PIR" = 6,
                                            "GO ID" = 7,
                                            "KEGG ID" = 8,
                                            "Pathway" = 9,
                                            "Unipathway" = 10),
                             selected = c(3:10), inline = TRUE)
          ),
        conditionalPanel(
          condition = "input.specie == 10090",
          checkboxGroupInput("mouse_var", label = "",
                             choices = list("UCSC" = 3,
                                            "PDB" = 4,
                                            "Reactome" = 5,
                                            "PIR" = 6,
                                            "GO ID" = 7,
                                            "KEGG ID" = 8,
                                            "Pathway" = 9,
                                            "Unipathway" = 10),
                             selected = c(3:10), inline = TRUE)
        ),
        conditionalPanel(
          condition = "input.specie == 3702",
          checkboxGroupInput("cress_var", label = "",
                             choices = list("PDB" = 3,
                                            "Reactome" = 4,
                                            "PIR" = 5,
                                            "GO ID" = 6,
                                            "KEGG ID" = 7,
                                            "Pathway" = 8,
                                            "Unipathway" = 9),
                             selected = c(3:9), inline = TRUE)
        ),
        conditionalPanel(
          condition = "input.specie == 7227",
          checkboxGroupInput("fly_var", label = "",
                             choices = list("UCSC" = 3,
                                            "PDB" = 4,
                                            "Reactome" = 5,
                                            "PIR" = 6,
                                            "GO ID" = 7,
                                            "KEGG ID" = 8,
                                            "Pathway" = 9,
                                            "Unipathway" = 10,
                                            "Flybase" = 11),
                             selected = c(3:11), inline = TRUE)
          ),
        conditionalPanel(
          condition = "input.specie == 6239",
          checkboxGroupInput("worm_var", label = "",
                             choices = list("UCSC" = 3,
                                            "PDB" = 4,
                                            "Reactome" = 5,
                                            "PIR" = 6,
                                            "GO ID" = 7,
                                            "KEGG ID" = 8,
                                            "Pathway" = 9,
                                            "Unipathway" = 10,
                                            "Wormbase" = 11),
                             selected = c(3:11), inline = TRUE)
          ),
        conditionalPanel(
          condition = "input.specie == 559292",
          checkboxGroupInput("yeast_var", label = "",
                             choices = list("PDB" = 3,
                                            "Reactome" = 4,
                                            "PIR" = 5,
                                            "GO ID" = 6,
                                            "KEGG ID" = 7,
                                            "Pathway" = 8,
                                            "Unipathway" = 9),
                             selected = c(3:9), inline = TRUE)
          )
        )
      ),
    hr(),
    mainPanel(
      width = 12,
      checkboxInput("sidebar", "Hide the panel", value = FALSE),
      tabsetPanel(
        tabPanel("Initial data",
                 h4("Initial table"),
                 dataTableOutput("Input")
        ),
        tabPanel("Statistics",
                 br(),
                 h4("Proportion test"),
                 textOutput("PValTest"),
                 h4("Conclusion"),
                 textOutput("Conclusion"),
                 br(),
                 br(),
                 br(),
                 textOutput("DimUpDown"),
                 downloadButton("downUpDownData", "here")
        ),
        tabPanel("Graphics",
                 fluidRow(
                   column(4, plotOutput("PieNew")),
                   column(4, plotOutput("Barplot")),
                   column(4, plotOutput("PieTotal"))
                 ),
                 fluidRow(
                   column(4, align = "center", downloadButton("PieNewDown", "Download graph")),
                   column(4, align = "center", downloadButton("BarplotDown", "Download graph")),
                   column(4, align = "center", downloadButton("PieTotalDown", "Download graph"))
                 ),
                 fluidRow(
                   column(8, plotOutput("PieDeg")),
                   column(4,
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          downloadButton("PieDegDown", "Download graph"))
                 )
                 
        ),
        tabPanel("Final results",
                 fluidRow(
                   column(4,
                          h4("Dimension"),
                          textOutput("DimFinalData")
                   ),
                   column(4,
                          br(),
                          downloadButton("downFinalData", "Download table")
                   )
                 ),
                 h4("Final data"),
                 dataTableOutput("FinalData")
        ),
        tabPanel("Differentially Expressed Genes analysis",
                   fluidRow(
                     column(6,
                            h4("Dimension"),
                            textOutput("DimRem_Use")),
                     column(6,
                            br(),
                            downloadButton("VolcanoDown", "Download Volcano plot")
                     )
                   ),
                   plotOutput("Volcano"),
                   fluidRow(
                     column(6, h4("DEG table")),
                     column(6,
                            br(),
                            downloadButton("downRem_Use", "Download table")
                     )
                   ),
                   dataTableOutput("Rem_Use")
        ),
        tabPanel("GO",
                 fluidRow(
                   column(6,
                          h4("Dimension"),
                          textOutput("DimDFGO")
                   ),
                   column(6,
                          br(),
                          downloadButton("downDFGO", "Download table")
                   )
                 ),
                 h4("GO table"),
                 dataTableOutput("DFGO")
        ),
        tabPanel("GO graphs",
                 br(),
                 downloadButton("downBarplotGO", "Download barplot"),
                 plotOutput("BarplotGO", height = "1000px")
        )
      )
    )
  )
)