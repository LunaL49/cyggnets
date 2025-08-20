library(shiny)
library(shinycssloaders)
library(shinyjs)
library(httr)
library(jsonlite)
source("compound2info.R")
source("plot_network.R")

file_names = readRDS("all_file_names.rds")
compounds = unique(sub("_.*", "", file_names))

base_download_url = "https://huggingface.co/datasets/lunaql/cyggnets/resolve/main/"

ui = fluidPage(
  useShinyjs(),
  
  # Custom CSS for full-height plots
  tags$head(
    tags$style(HTML("
      body {
        margin: 0;
        padding: 0;
      }
      
      .container-fluid {
        padding-left: 0 !important;
        padding-right: 0 !important;
        margin-left: 12px !important;
        margin-right: 0 !important;
        max-width: 100% !important;
        width: 100% !important;
      }
      
      .row {
        margin-left: 0 !important;
        margin-right: 0 !important;
      }
      
      .info-box {
        max-width: 100%;
        word-wrap: break-word;
        margin-bottom: 10px;
      }
      
      .full-height-plots {
        height: calc(100vh - 200px);
      }
      
      .full-height-plots .shiny-output-container {
        height: 100%;
        width: 100%;
      }
      
      #left_plot, #right_plot {
        height: calc(100vh - 200px) !important;
      }
    "))
  ),
  
  # Header Title
  titlePanel("Cyggnets - CycleGAN RNAseq for drug and target network signatures"),
  
  h5(HTML('Cyggnets is trying to answer the question: how similar are the transcriptomic effects of a chemical perturbation
     compared to a biological perturbation at the same target gene/protein? <br>
     The data is based on LINCS L1000, which has been synthetically amplified to full RNAseq by Avi Ma\'ayan\'s group in 
     their 2022 paper "Transforming L1000 profiles to RNA-seq-like profiles with deep learning" (Jeon et al.). <br>
     All precomputed GSEA files are hosted on Hugging Face, publicly accessible at: <a href="https://huggingface.co/datasets/lunaql/cyggnets/tree/main">
     https://huggingface.co/datasets/lunaql/cyggnets/tree/main</a>.')),
  
  # Layout with sidebar and main panel
  sidebarLayout(
    
    # Left panel
    sidebarPanel(
      width = 3,
      
      selectInput(inputId = "compound_name", label = "First, choose a compound to start:", 
                  choices = c("Choose a compound..." = "", compounds),
                  selected = ""),
      
      withSpinner(uiOutput("info")),
      
      uiOutput("ui_input2"),
      
      uiOutput("ui_input3"),
      
      uiOutput("ui_input4"),
      
      uiOutput("ui_input5"),
      
      actionButton("compute", "Compute")
    ),
    
    # Right panel (for plots or outputs)
    mainPanel(
      width = 9,
      
      selectInput("plot_choice", "Select plot pair to view:",
                  choices = c("Hallmark gene sets", "Canonical pathways gene sets (top 50)"),
                  selected = "Hallmark gene sets"),
      
      p("Left panel - chemical perturbation; right panel - biological perturbation. Red - upregulated; 
        blue - downregulated. Feel free to click, pan, zoom in, or use the 'select by ID' dropdown 
        to highlight pathways of interest."),
      
      # Plots with full height
      div(class = "full-height-plots",
          fluidRow(
            column(6, visNetworkOutput("left_plot")),
            column(6, visNetworkOutput("right_plot"))
          )
      )
    )
  )
)

server = function(input, output, session) {
  
  output$info = renderUI({
    req(input$compound_name != "")
    response = compound_info(input$compound_name)
    indications = response$indications[1:min(5, length(response$indications))]
    indications = tolower(indications)
    mechanisms = response$mechanisms[1:min(3, length(response$mechanisms))]
    
    # Visible box with styling
    tagList(
      tags$p(paste("Diseases where",input$compound_name, "is used in includes:")),
      
      tags$div(class = "info-box",
               style = "border: 1px solid #ccc; padding: 10px; background-color: #f9f9f9;",
               HTML(paste(indications, collapse = ", "))
      ),
      
      tags$p(paste("It has these known mechanism-of-actions:")),
      
      tags$div(class = "info-box",
               style = "border: 1px solid #ccc; padding: 10px; background-color: #f9f9f9;",
               HTML(paste(mechanisms, collapse = ", "))
      ),
      
      tags$p("These are its known targets in the human genome:"),
      
      tags$div(class = "info-box",
               style = "border: 1px solid #ccc; padding: 10px; background-color: #f9f9f9;",
               HTML(paste(response$target_genes, collapse = ", "))
      )
    )
  })
  
  file_names_1 = reactive({
    file_names[grepl(paste0("^", input$compound_name), file_names)]
  })
  
  # Render second dropdown
  output$ui_input2 = renderUI({
    req(input$compound_name != "")
    req(file_names_1()) 
    targets = sub("^[^_]*_", "", file_names_1())
    targets = unique(sub("_.*", "", targets))
    selectInput(inputId = "target_name", 
                label = "Choose an available target gene:",
                choices = c("Choose a target gene..." = "", targets),
                selected = "")
  })
  
  file_names_2 = reactive({
    req(file_names_1())
    req(input$target_name != "")
    file_names_2 = sub("^[^_]*_", "", file_names_1())
    file_names_2 = file_names_2[grepl(paste0("^", input$target_name), file_names_2)]
  })
  
  output$ui_input3 = renderUI({
    req(input$target_name != "")
    req(file_names_2()) 
    perts = sub("^[^_]*_", "", file_names_2())
    perts = unique(sub("_.*", "", perts))
    selectInput(inputId = "pert_name", 
                label = "Choose an available biological perturbation method:",
                choices = c("Choose a perturbation method..." = "", perts),
                selected = "")
  })
  
  file_names_3 = reactive({
    req(file_names_2())
    req(input$pert_name != "")
    file_names_3 = sub("^[^_]*_", "", file_names_2())
    file_names_3 = file_names_3[grepl(paste0("^", input$pert_name), file_names_3)]
  })
  
  output$ui_input4 = renderUI({
    req(input$pert_name != "")
    req(file_names_3()) 
    cells = sub("^[^_]*_", "", file_names_3())
    cells = unique(sub("_.*", "", cells))
    selectInput(inputId = "cell_name", 
                label = "Choose an available cell line:",
                choices = c("Choose a cell line..." = "", cells),
                selected = "")
  })
  
  file_names_4 = reactive({
    req(file_names_3())
    req(input$cell_name != "")
    file_names_4 = sub("^[^_]*_", "", file_names_3())
    file_names_4 = file_names_4[grepl(paste0("^", input$cell_name), file_names_4)]
  })
  
  output$ui_input5 = renderUI({
    req(input$cell_name != "")
    req(file_names_4()) 
    doses = sub("^[^_]*_", "", file_names_4())
    doses = sub("\\.rds$", "", doses)
    doses = unique(sub("_.*", "", doses))
    selectInput(inputId = "dose_name", 
                label = "Choose a compound dose (if available, try 10 uM first):",
                choices = c("Choose a dose..." = "", doses),
                selected = "")
  })
  
  all_options_selected = reactive({
    all(
      !is.null(input$compound_name) && input$compound_name != "",
      !is.null(input$target_name) && input$target_name != "",
      !is.null(input$pert_name) && input$pert_name != "",
      !is.null(input$cell_name) && input$cell_name != "",
      !is.null(input$dose_name) && input$dose_name != ""
    )
  })
  
  observe({
    if (all_options_selected()) {
      enable("compute")
    } else {
      disable("compute")
    }
  })
  
  observe({
    toggleState("plot_choice", condition = input$compute > 0)
  })
  
  data_to_plot = eventReactive(input$compute, {
    req(input$compound_name, input$target_name, input$pert_name, 
        input$cell_name, input$dose_name)
    filename = paste(input$compound_name, input$target_name, 
                     input$pert_name, input$cell_name, input$dose_name, sep = "_")
    filename = paste0(filename, ".rds")
    filename = gsub(" ", "%20", filename)
    file_download_url = paste0(base_download_url, filename, "?download=true")
    con = url(file_download_url)
    downloaded_data = readRDS(con)
    close(con)
    
    return(downloaded_data)
  })
  
  output$left_plot = renderVisNetwork({
    req(input$compute)
    data = data_to_plot()
    
    data_left = switch(input$plot_choice,
                       "Hallmark gene sets" = data$cp$GSEA_H,
                       "Canonical pathways gene sets (top 50)" = data$cp$GSEA_C2CP)
    
    plot_network(data_left)
  })
  
  output$right_plot = renderVisNetwork({
    req(input$compute)
    data = data_to_plot()
    
    data_right = switch(input$plot_choice,
                        "Hallmark gene sets" = data$bp$GSEA_H,
                        "Canonical pathways gene sets (top 50)" = data$bp$GSEA_C2CP)
    
    plot_network(data_right)
  })
  
}

shinyApp(ui, server)