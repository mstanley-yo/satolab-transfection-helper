# satolab-transfection-helper
# Load R packages
library(shiny)
library(bslib)
library(shinythemes)
library(tidyverse)
library(flextable)
library(officer)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "Transfection calculator",
    tabPanel(
      "Calculator",
      sidebarPanel(
        # input names and concentrations
        textInput("sample_name_input", "Sample name", value = NA),
        numericInput("spike_input", "Spike concentration (ng/µl)", value = NA, min = 0), # test: 500
        numericInput("hibit_input", "HiBiT concentration (ng/µl)", value = NA, min = 0), # test: 953
        numericInput("luc2_input", "Luc2 concentration (ng/µl)", value = NA, min = 0), # test: 920

        # select what type of well to transfect, and how many.
        layout_columns(
          radioButtons(
            "plate_input",
            "Plate/Dish type",
            list("6-well (2 mL)" = 2, "12-well (1 mL)" = 1, "15 cm (20 mL)" = 20)
          ),
          numericInput("num_to_transfect_input", "Number to transfect", value = NA, min = 0) # test: 3
        ),
        br(),

        # submit sample
        actionButton("add_sample", "Add Sample", icon = icon("plus")),
        br(), br(),
        actionButton("remove_all_samples", "Remove All Samples", icon = icon("trash"))
      ), # sidebarPanel
      mainPanel(
        uiOutput("pretty_output"),
        br(),
        p(HTML("INFO:<br>
                 Mass of DNA required for transfection is 1 µg per mL of medium.<br>
                 For pseudoviruses, the DNA mix is 40% HiBiT, 40% Luc2, and 20% S.")),
        br(),
        downloadButton("download_docx", "Download table as .docx")
      ) # mainPanel
    ), # Navbar, tabPanel Calculator
    tabPanel(
      "About",
      mainPanel(
        h3("About this website"),
        p("Written in R Shiny by Maximilian Stanley Yo."),
        p(
          "Follow development here: ",
          tags$a("GitHub Repository", href = "https://github.com/mstanley-yo/satolab-transfection-helper", target = "_blank")
        )
      )
    ) # Navbar, tabPanel About
  ) # navbarPage
) # fluidPage


# Define server function
server <- function(input, output) {
  # Function to validate inputs
  validateInputs <- function() {
    validate(
      need(input$sample_name_input != "", "Please enter a sample name."),
      need(!is.na(input$spike_input), "Please enter Spike concentration."),
      need(!is.na(input$hibit_input), "Please enter HiBiT concentration."),
      need(!is.na(input$luc2_input), "Please enter Luc2 concentration."),
      need(!is.na(input$num_to_transfect_input), "Please enter the number of plates to transfect.")
    )
  }

  # Set up dynamic sample table
  sample_data <- reactiveVal(data.frame(
    sample_id = character(),
    conc_spike = numeric(),
    volume_spike = numeric(),
    volume_hibit = numeric(),
    volume_luc2 = numeric(),
    volume_optimem = numeric(),
    volume_transit = numeric(),
    volume_transfect = numeric() 
  ))

  # add sample
  observeEvent(input$add_sample, {
    # validate that all inputs are present
    validateInputs()

    # Calculate total cell medium volume based on plate type and number of wells/dishes
    volume_cell_medium <- as.numeric(input$plate_input) * as.numeric(input$num_to_transfect_input)

    # Calculate reagent volumes to add (µL) based on desired mass and user-provided concentrations
    volume_spike_calc <- (200 * volume_cell_medium) / as.numeric(input$spike_input) # Spike: desired mass / concentration
    volume_hibit_calc <- (400 * volume_cell_medium) / as.numeric(input$hibit_input) # HiBiT
    volume_luc2_calc <- (400 * volume_cell_medium) / as.numeric(input$luc2_input) # Luc2

    # Transfection reagents (converted from mL into uL)
    volume_optimem_calc <- volume_cell_medium * 100
    volume_transit_calc <- volume_cell_medium * 3

    # consolidate into new row
    new_row <- data.frame(
      sample_id = input$sample_name_input,
      conc_spike = input$spike_input,
      volume_spike = volume_spike_calc,
      volume_hibit = volume_hibit_calc,
      volume_luc2 = volume_luc2_calc,
      volume_optimem = volume_optimem_calc,
      volume_transit = volume_transit_calc,
      volume_transfect = volume_cell_medium
    )

    # bind to current
    current <- sample_data()
    sample_data(rbind(current, new_row))
  }) # observeEvent - add_sample

  # remove all samples
  observeEvent(input$remove_all_samples, {
    # set sample_data to empty.
    sample_data(data.frame(
      sample_id = character(),
      conc_spike = numeric(),
      volume_spike = numeric(),
      volume_hibit = numeric(),
      volume_luc2 = numeric(),
      volume_optimem = numeric(),
      volume_transit = numeric(),
      volume_transfect = numeric()
    ))
  }) # observeEvent - remove_all_samples

  # reactive df for display and for creating flextable. Rename column headers and round values
  df_display <- reactive({
    sample_data() %>%
      mutate(across(where(is.numeric), function(x) round(x, 2))) %>%
      rename(
        Sample = sample_id,
        `S conc. (ng/µL)` = conc_spike,
        `S volume (µL)` = volume_spike,
        `HiBiT volume (µL)` = volume_hibit,
        `Luc2 volume (µL)` = volume_luc2,
        `Add OptiMEM (µL) ` = volume_optimem,
        `Add TransIT (µL)` = volume_transit,
        `Transfect to: (mL)` = volume_transfect
      )
  })

  output$sample_table <- renderTable({
    validateInputs() # validate that all inputs are present
    df_display() # render table
  }) # renderTable

  # flextable reactive object.
  ft <- reactive({
    df_display() %>%
      flextable() %>%
      autofit() %>%
      add_header_lines(
        values = paste0("Pseudovirus transfection (", Sys.Date(), ") dilution table")
      ) %>%
      add_footer_lines(
        values = as_paragraph(
          paste0("HiBiT plasmid concentration: ", input$hibit_input, "ng/µL\n"),
          paste0("Luc2 plasmid concentration: ", input$luc2_input, "ng/µL\n"),
          "\n",
          "Protocol:\n",
          "1. Add calculated volumes of DNA to a sterile tube.\n",
          "2. Add calculated volume of Opti-Mem.\n",
          "3. Add calculated volume of Trans-IT.\n",
          "4. Incubate samples for 15 minutes at room temperature.\n",
          "5. Add samples to cell medium at 10% v/v. For example, add 200 µL to each 2 mL well."
        )
      ) %>%
      bold(bold = TRUE, part = "header") %>%
      bold(j = 3:5, bold = TRUE, part = "body") %>%
      color(j = 3:5, color = "red", part = "body") %>%
      align(align = "center", part = "header") %>%
      align(align = "center", part = "body")
  })

  output$pretty_output <- renderUI({
    # validate that all inputs are present
    validateInputs()

    # Render flextable as HTML widget
    tagList(
      h3(paste0("Transfection protocol (", Sys.Date(), ")")),
      flextable::htmltools_value(ft())
    )
  })

  # render flextable as .docx and make downloadable
  output$download_docx <- downloadHandler(
    filename = function() {
      paste0("transfection_table_", Sys.Date(), ".docx")
    },
    content = function(file) {
      ft_width <- 9 # increase to extend the width of the flextable
      ft_docx <- ft() %>%
        width(width = dim(.)$widths * ft_width / (flextable_dim(.)$widths)) # format for docx.

      # Save as docx
      sect_properties <- prop_section(
        page_size = page_size(
          orient = "landscape",
          width = 7, height = 10
        ),
        type = "continuous",
        page_margins = page_mar(
          top = 0.5, right = 0.5, bottom = 0.5, left = 0.5, header = 0.3, footer = 0.3, gutter = 0
        )
      )
        
      flextable::save_as_docx(ft_docx, path = file, pr_section = sect_properties)
    }
  )
} # server

# Create Shiny object
shinyApp(ui = ui, server = server)
