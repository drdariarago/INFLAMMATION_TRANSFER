# shiny test

# Load libraries 

library(shiny)
library(tidyverse)
library(here)

## Load dataset
full_data <- 
  read_csv("../results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>%
  mutate(
    timepoint = factor( x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  ) %>%
  select(-contrast) %>%
  pivot_wider(names_from = exposure, values_from = c(q_value, logFC)) %>%
  rename(log_base_counts = logFC_baseline, log_fc_response = logFC_response)

### Define page layout
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("regexp",
                 textInput("pattern", label = "Gene ID (as regular expression)", placeholder = "Adam1"),
                 actionButton(inputId = "run_regex", label = "Search")
        ),
        tabPanel("genelist",
                 textAreaInput("genelist", label = "list of gene IDs", placeholder = "Adam10\nAdam11"),
                 actionButton(inputId = "run_genelist", label = "Search")
        )
      ),
      tags$hr(),
      numericInput("q_val", label = "q value threshold", value = 0.05),
      checkboxGroupInput("tissues", "Choose tissues:",
                         choiceNames = list("lung", "maternal liver", "placentas", "foetal liver"),
                         choiceValues = list("maternal_lung", "maternal_liver", "placentas", "fetal_liver")),
      tags$hr(),
      tags$h4("Genes found"),
      textOutput(outputId = "genelist")
    ),
    mainPanel(
      plotOutput("time_series")
      # plotOutput("ma")
    )
  )
)

server <- function(input, output){

## Filter genes of interest

  rna_data <- reactiveValues(data = NULL)
  
  observeEvent( input$run_regex, 
                rna_data$data <- 
                  full_data %>% 
                  filter( tissue %in% input$tissues ) %>%
                  filter( grepl(isolate(input$pattern), mgi_symbol) )
                )
  
  observeEvent( input$run_genelist,
                rna_data$data <-
                  full_data %>% 
                  filter(tissue %in% input$tissues ) %>% 
                  filter( 
                    mgi_symbol %in% (input$genelist %>% strsplit(., "\\s+") %>% extract2(1) )
                  )
  )

## List of genes found
output$genelist <- renderText(
  rna_data$data %>%
    pull(mgi_symbol) %>%
    unique()
)
  
## Plot results
  output$time_series <- renderPlot(
    rna_data$data %>%
      ggplot(
        aes(
          x = timepoint, y = log_fc_response,
          shape = q_value_response < isolate(input$q_val),
          group = mgi_symbol, label = mgi_symbol
        )
      ) +
      geom_line(col = 'black', alpha = 0.5) +
      geom_point(aes(col = log_base_counts)) +
      facet_wrap( ~ tissue) +
      viridis::scale_color_viridis() +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      labs(shape = "significant", col = "log base counts") +
      ggrepel::geom_text_repel(
        data = rna_data$data %>%
          filter(q_value_response < isolate(input$q_val) ),
        aes(
          x = timepoint, y = log_fc_response, label = mgi_symbol
        )
      ) +
      ggtitle("Log fold change over time")
  )
  
## MA plots
  output$ma <- renderPlot(
    rna_data$data %>% 
      ggplot(
        aes(x = log_base_counts, y = log_fc_response,
            col = q_value_response < isolate(input$q_val), 
            label = mgi_symbol)
      ) +
      geom_point() +
      ggrepel::geom_text_repel(
        data = rna_data$data %>%
          filter(q_value_response < isolate(input$q_val) ),
        aes(
          x = log_base_counts, y = log_fc_response, 
          label = mgi_symbol
        )
      ) +
      facet_grid(timepoint ~ tissue) +
      labs(color = "significant") +
      theme_minimal() +
      ggtitle("MA plot") 
  )
}

shinyApp(ui = ui, server = server)
