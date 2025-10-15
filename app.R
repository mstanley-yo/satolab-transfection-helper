# satolab-transfection-helper
# Load R packages
library(shiny)
library(bslib)
library(shinythemes)
library(tidyverse)
library(flextable)
library(officer)

# define variables
github_link <- "https://github.com/mstanley-yo/satolab-transfection-helper"

# function to round volumes based on if using p20, p200, or p1000
round_volumes <- function(vol) {
    case_when(
        vol < 20 ~ round(vol, 2),
        vol > 200 ~ round(vol),
        .default = round(vol, 1)
    )
}

# ui function #####
ui <- page_fluid(
    # theme
    theme = bs_theme(bootswatch = "flatly"),
    
    # title
    tags$h3(
        "Pseudovirus Transfection Calculator", 
        class = "text-primary", 
        style = "margin-top: 15px;margin-bottom: 15px;"
    ),
    
    # sidebar layout
    sidebarLayout(
        sidebarPanel(
            textInput(
                "sample_name_input",
                "Sample name", 
                value = NA
            ),
            numericInput(
                "spike_input", 
                "Spike concentration (ng/µL)", 
                value = NA, 
                min = 0
            ),
            numericInput(
                "hibit_input", 
                "HiBiT concentration (ng/µL)", 
                value = NA, 
                min = 0
            ),
            numericInput(
                "luc2_input", 
                "Luc2 concentration (ng/µL)", 
                value = NA, 
                min = 0
            ),
            
            layout_columns(
                radioButtons(
                    "plate_input",
                    "Plate/Dish type",
                    choices = list(
                        "6-well (2 mL)" = 2, 
                        "12-well (1 mL)" = 1, 
                        "15 cm (20 mL)" = 20
                    ),
                    selected = 20
                ),
                numericInput(
                    "num_input", 
                    "Number to transfect", 
                    value = NA, 
                    min = 0
                )
            ),
            
            br(),
            actionButton(
                "add_sample", 
                "Add Sample", 
                icon = icon("plus"),
                class = "btn-primary"
            ),
            br(), 
            br(),
            actionButton(
                "remove_all_samples", 
                "Remove All Samples", 
                icon = icon("trash"),
                class = "btn-danger"
            )
        ),
        
        mainPanel(
            uiOutput("pretty_output"),
            br(),
            downloadButton("download_docx", "Download table as .docx"),
            br(),
            br(),
            p("Written in R Shiny by Maximilian Stanley Yo."),
            p(
                "Follow development here: ",
                tags$a(
                    "GitHub Repository", 
                    href = github_link,
                    target = "_blank"
                )
            )
        )
    )
)

# server function #####
server <- function(input, output) {

    validate_inputs <- function() {
        # display anyway if sample table is not empty
        if (nrow(sample_data()) > 0) {
                return(invisible(TRUE))
        }
        
        # else validate inputs
        validate(
            need(
                input$sample_name_input != "", 
                "Please enter a sample name."
            ),
            need(
                !is.na(input$spike_input), 
                "Please enter Spike concentration."
            ),
            need(
                !is.na(input$hibit_input), 
                "Please enter HiBiT concentration."),
            need(
                !is.na(input$luc2_input), 
                "Please enter Luc2 concentration."
            ),
            need(
                !is.na(input$num_input), 
                "Please enter the number of plates to transfect."
            )
        )
    }
    
    # setup reactive sample table
    sample_data <- reactiveVal(
        data.frame(
            sample_id = character(),
            conc_spike = numeric(),
            volume_spike = numeric(),
            volume_hibit = numeric(),
            volume_luc2 = numeric(),
            volume_master = numeric(),
            volume_optimem = numeric(),
            volume_transit = numeric(),
            volume_transfect = numeric() 
        )
    )
    
    # display reactive sample table as flextable
    ft <- reactive({
        # calculate volumes
        vols <- calc_volumes()
        
        # rename column headers and round values
        data <- sample_data() %>%
            mutate(across(where(is.numeric), round_volumes)) %>%
            rename(
                Sample = sample_id,
                `S conc. (ng/µL)` = conc_spike,
                `S volume (µL)` = volume_spike,
                `HiBiT volume (µL)` = volume_hibit,
                `Luc2 volume (µL)` = volume_luc2,
                `Master mix volume (uL)` = volume_master,
                `Add OptiMEM (µL) ` = volume_optimem,
                `Add TransIT (µL)` = volume_transit,
                `Transfect to: (mL)` = volume_transfect
            )
        
        # master mix text
        mastermix_text <- paste(
            "To create a master mix, add",
            round_volumes(
                vols$volume_hibit_calc * (nrow(sample_data()) + 0.5)
            ),
            "uL HiBiT plasmid and",
            round_volumes(
                vols$volume_luc2_calc * (nrow(sample_data()) + 0.5)
            ),
            "uL Luc2 plasmid to a master mix tube,",
            "then add the indicated master mix volume to each tube.\n"
        )
        
        # protocol text
        protocol_text <- paste(
            "Protocol:\n",
            "1. Add calculated volumes of DNA to a sterile tube.\n",
            "2. Add calculated volume of Opti-Mem.\n",
            "3. Add calculated volume of Trans-IT.\n",
            "4. Incubate samples for 15 minutes at room temperature.\n",
            "5. Add sample volume equivalent to 10% of cell medium volume. 
            (For example, add 2 mL to each 20 mL plate.)"
        )
        
        # process into flextable
        data %>%
            flextable() %>%
            autofit() %>%
            add_header_lines(
                values = paste0(
                    "Pseudovirus transfection (",
                    Sys.Date(), 
                    ") dilution table"
                )
            ) %>%
            add_footer_lines(
                values = as_paragraph(
                    paste(
                        "HiBiT plasmid concentration:", 
                        input$hibit_input, 
                        "ng/µL\n"
                    ),
                    paste(
                        "Luc2 plasmid concentration:", 
                        input$luc2_input, 
                        "ng/µL\n"
                    ),
                    "\n",
                    mastermix_text,
                    "\n",
                    protocol_text
                )
            ) %>%
            bold(bold = TRUE, part = "header") %>%
            bold(j = 3:6, bold = TRUE, part = "body") %>%
            color(j = 3:5, color = "red", part = "body") %>%
            color(j = 6, color = "blue", part = "body") %>%
            align(align = "center", part = "header") %>%
            align(align = "center", part = "body")
    })
    
    # Reactive expression for calculated volumes
    calc_volumes <- reactive({
        validate_inputs()
        
        # Calculate total cell medium volume
        volume_cell_medium <- as.numeric(input$plate_input) * input$num_input
        
        # Calculate reagent volumes
        volume_spike_calc <- (200 * volume_cell_medium) / input$spike_input
        volume_hibit_calc <- (400 * volume_cell_medium) / input$hibit_input
        volume_luc2_calc <- (400 * volume_cell_medium) / input$luc2_input
        volume_optimem_calc <- volume_cell_medium * 100
        volume_transit_calc <- volume_cell_medium * 3
        volume_master_calc <- volume_hibit_calc + volume_luc2_calc
        
        # Return as named list
        list(
            volume_cell_medium = volume_cell_medium,
            volume_spike_calc = volume_spike_calc,
            volume_hibit_calc = volume_hibit_calc,
            volume_luc2_calc = volume_luc2_calc,
            volume_master_calc = volume_master_calc,
            volume_optimem_calc = volume_optimem_calc,
            volume_transit_calc = volume_transit_calc
        )
    })
    
    # button - add sample
    observeEvent(input$add_sample, {
            validate_inputs()
            
            vols <- calc_volumes()
            
            # consolidate into new row
            new_row <- data.frame(
                    sample_id = input$sample_name_input,
                    conc_spike = input$spike_input,
                    volume_spike = vols$volume_spike_calc,
                    volume_hibit = vols$volume_hibit_calc,
                    volume_luc2 = vols$volume_luc2_calc,
                    volume_master = vols$volume_master_calc,
                    volume_optimem = vols$volume_optimem_calc,
                    volume_transit = vols$volume_transit_calc,
                    volume_transfect = vols$volume_cell_medium
            )
            
            current <- sample_data()
            sample_data(rbind(current, new_row))
    }) # observeEvent - add_sample

    # button - remove all samples
    observeEvent(input$remove_all_samples, {
        # set sample_data to empty.
        sample_data(data.frame(
            sample_id = character(),
            conc_spike = numeric(),
            volume_spike = numeric(),
            volume_hibit = numeric(),
            volume_luc2 = numeric(),
            volume_master = numeric(),
            volume_optimem = numeric(),
            volume_transit = numeric(),
            volume_transfect = numeric()
        ))
    }) # observeEvent - remove_all_samples
    
    # render flextable as HTML widget
    output$pretty_output <- renderUI({
        validate_inputs()
        
        tagList(
            h4(paste0("Transfection protocol (", Sys.Date(), ")")),
            flextable::htmltools_value(ft())
        )
    })

    # render flextable as .docx and make downloadable
    output$download_docx <- downloadHandler(
        filename = function() {
            paste0("transfection_table_", Sys.Date(), ".docx")
        },
        content = function(file) {
            # format for docx.
            ft_width <- 9 # increase to extend the width of the flextable
            ft_docx <- ft() %>%
                width(width = dim(.)$widths * ft_width /
                          (flextable_dim(.)$widths)) 

            # save as docx
            sect_properties <- prop_section(
                page_size = page_size(
                    orient = "landscape",
                    width = 7, height = 10
                ),
                type = "continuous",
                
                # make margins narrow
                page_margins = page_mar(
                    top = 0.5, 
                    right = 0.5, 
                    bottom = 0.5, 
                    left = 0.5, 
                    header = 0.3, 
                    footer = 0.3, 
                    gutter = 0
                )
            )
            flextable::save_as_docx(
                    ft_docx, 
                    path = file, 
                    pr_section = sect_properties
            )
        }
    )
}


# Create Shiny object #####
shinyApp(ui = ui, server = server)
