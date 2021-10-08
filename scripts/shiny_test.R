# shiny test

# Load libraries 

library(shiny)
library(tidyverse)
library(magrittr)
library(here)

## Load datasets
# Limma results
limma_data <- 
  read_csv("../results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>%
  mutate(
    timepoint = factor( x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  ) %>%
  select(-contrast) %>%
  pivot_wider(names_from = exposure, values_from = c(q_value, logFC)) %>%
  rename(log_base_counts = logFC_baseline, log_fc_response = logFC_response)

# Raw TPM
tpm_data <-
  here::here("results/tpm_summary/tpm_tibble.rds") %>% 
  read_rds()

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
      tabsetPanel(
        tabPanel("Fold Change",
                 plotOutput("time_series")
        ),
        tabPanel("TPM time series",
                 plotOutput("tpm_plot")
        ),
        tabPanel("TPM heatmap",
                 plotOutput("heatmap"))
      )
    )
  )
)

server <- function(input, output){

## Filter genes of interest

  rna_data <- reactiveValues(logfc = NULL)
  
  observeEvent( input$run_regex, 
                rna_data$logfc <- 
                  limma_data %>% 
                  filter( tissue %in% input$tissues ) %>%
                  filter( grepl(isolate(input$pattern), mgi_symbol) )
                )

  observeEvent( input$run_genelist,
                rna_data$logfc <-
                  limma_data %>% 
                  filter(tissue %in% input$tissues ) %>% 
                  filter( 
                    mgi_symbol %in% ( 
                      input$genelist %>% 
                        strsplit(., "\\s+") %>% 
                        unlist
                    )
                  )
  )

## List of genes found
output$genelist <- renderText(
  rna_data$logfc %>%
    pull(mgi_symbol) %>%
    unique()
)
  
## Fold change over time
  output$time_series <- renderPlot(
    rna_data$logfc %>%
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
        data = rna_data$logfc %>%
          filter(q_value_response < isolate(input$q_val) ),
        aes(
          x = timepoint, y = log_fc_response, label = mgi_symbol
        )
      ) +
      ggtitle("Log fold change over time")
  )
  
## TPM over time
  
  tpm_plot_data <- reactiveValues(tpm = NULL)
  
  observeEvent( input$run_regex,
                tpm_plot_data$tpm <-
                  tpm_data %>%
                  filter( tissue %in% input$tissues ) %>%
                  filter( grepl(isolate(input$pattern), mgi_symbol) )
  )
  
  observeEvent( input$run_genelist,
                tpm_plot_data$tpm <-
                  tpm_data %>% 
                  filter(tissue %in% input$tissues ) %>% 
                  filter( 
                    mgi_symbol %in% (input$genelist %>% strsplit(., "\\s+") %>% extract2(1) )
                  )
  )
  
  output$tpm_plot <- renderPlot(
    tpm_plot_data$tpm %>%
      ggplot(
        aes(
          x = as.factor(timepoint), y = log10(TPM+1),
          group = exposure, col = exposure
        )
      ) +
      geom_boxplot(aes(group = timepoint)) +
      geom_jitter(height = 0, alpha = 0.7) +
      facet_grid(tissue ~ mgi_symbol, scales = "free") +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      ggtitle("TPM over time") 
  )
  
  output$heatmap <- renderPlot(
    tpm_plot_data$tpm %>% 
      select(sample_id, mgi_symbol, TPM) %>% 
      pivot_wider(
        id_cols = sample_id, 
        names_from = mgi_symbol, 
        values_from = TPM
      ) %>% 
      column_to_rownames("sample_id") %>% 
      as.matrix() %>% 
      t() %>% 
      add(1) %>% 
      log10() %>% 
      pheatmap::pheatmap(
      mat = ., scale = "none",
      annotation_col = tpm_plot_data$tpm %>%
        select(sample_id, timepoint, exposure, tissue) %>% 
        distinct() %>% 
        column_to_rownames("sample_id"),
      show_colnames = FALSE,
      annotation_colors = list(
        exposure = c(LPS = "orange", ctr = "blue")
      )
    )
  )
  
}

shinyApp(ui = ui, server = server)
