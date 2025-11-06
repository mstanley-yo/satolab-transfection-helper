library(shiny)
library(bslib)
library(shinyvalidate)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)

# clickable github link + icon
github_link <- tags$a(
    shiny::icon("github"), "GitHub",
    href = "https://github.com/mstanley-yo/satolab-transfection-helper",
    target = "_blank"
)

# ui function #####
ui <- page_navbar(
    theme = bs_theme(bootswatch = "flatly"),
    title = "Pseudovirus Transfection Calculator",
    
    # calculator panel
    nav_panel(
        "Calculator",
        layout_columns(
            col_widths = c(4, 8),
            card(
                card_header("Add Transfection Sample"),
                layout_columns(
                    col_widths = c(12, 12),
                    textInput(
                        "sample_name",
                        "Spike", 
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
                    
                ),
                layout_columns(
                    radioButtons(
                        "plate_input",
                        "Plate/Dish type",
                        choices = list(
                            "6-well (per well, 2 mL)" = 2,
                            "12-well (per well, 1 mL)" = 1,
                            "15 cm dish (20 mL)" = 20
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
                actionButton(
                    "add_sample", 
                    "Add Sample", 
                    icon = icon("plus"),
                    class = "btn-primary"
                ),
                actionButton(
                    "remove_all_samples", 
                    "Remove All Samples", 
                    icon = icon("trash"),
                    class = "btn-danger"
                )
            ),
            
            card(
                card_header("Transfection Table & Protocol"),
                uiOutput("table_output"),
                downloadButton("download_docx", "Download table as .docx"),
                p("Written in R Shiny by Maximilian Stanley Yo."),
                github_link
            )
        )
    ),
    
    # settings panel
    nav_panel(
        "Settings",
        h4("Master mix settings"),
        checkboxInput("show_hibit_luc2", "Display HiBiT/Luc2 volumes", TRUE),
        checkboxInput("show_mm", "Combine HiBiT/Luc2 as master mix", TRUE),
        numericInput(
            "mm_extra", 
            "Extra HiBiT/Luc2 master mix (plates)", 
            "0.5", 
            min = 0
        ),
    )
)

# server function #####
server <- function(input, output) {
    validate_inputs <- function() {
        # display anyway if sample table is not empty
        if (nrow(sample_data()) > 0) {
                return(invisible(TRUE))
        }
        
        # else validate and require inputs
        validate(
            need(input$sample_name, "Please enter sample name."),
            need(input$spike_input, "Please enter Spike concentration."),
            need(input$hibit_input, "Please enter HiBiT concentration."),
            need(input$luc2_input, "Please enter Luc2 concentration."),
            need(input$num_input, "Please enter number of plates.")
        )
    }
    
    # validate inputs through warnings
    iv <- InputValidator$new()
    iv$add_rule("mm_extra", sv_required())
    iv$add_rule("mm_extra", function(value) {
        if (value < 0) {
            "Extra master mix (plates) should be at least 0."
        }
    })
    iv$add_rule("num_input", sv_optional())
    iv$add_rule("num_input", function(value) {
        if (value > 20) {
            "Volume may exceed 50 mL tube capacity. Should keep under 21."
        }
    })
    iv$enable()
    
    # setup reactive sample table
    empty_table <- data.frame(
        sample_id = character(),
        conc_spike = numeric(),
        volume_spike = numeric(),
        volume_hibit = numeric(),
        volume_luc2 = numeric(),
        volume_master = numeric(),
        volume_optimem = numeric(),
        volume_transit = numeric(),
        volume_transfect = numeric(),
        plate_count = character()
    )
    
    sample_data <- reactiveVal(
        empty_table
    )
    
    # function to round volumes based on if using p20, p200, or p1000
    round_volumes <- function(vol) {
        case_when(
            vol < 20 ~ round(vol, 2),
            vol > 200 ~ round(vol),
            .default = round(vol, 1)
        )
    }
    
    # Reactive expression for calculated volumes
    calc_volumes <- reactive({
        validate_inputs()
        
        # Calculate total cell medium volume
        volume_cell_medium <- as.numeric(input$plate_input) * input$num_input
        plate_count <- paste(
            input$num_input,
            case_when(
                input$plate_input == 2 ~ "6-well wells",
                input$plate_input == 1 ~ "12-well wells", 
                input$plate_input == 20 ~ "15-cm dishes"
            )
        )
        
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
            volume_transit_calc = volume_transit_calc,
            plate_count = plate_count
        )
    })
    
    # display reactive sample table as flextable
    ft <- reactive({
        # calculate volumes
        vols <- calc_volumes()
        
        # rename column headers and round values
        data <- sample_data() %>%
            rename(
                S = sample_id,
                `S conc.\n(ng/µL)` = conc_spike,
                `S vol.\n(µL)` = volume_spike,
                `HiBiT\n(µL)` = volume_hibit,
                `Luc2\n(µL)` = volume_luc2,
                `Master mix\n(µL)` = volume_master,
                `Opti-MEM\n(µL) ` = volume_optimem,
                `TransIT\n(µL)` = volume_transit,
                `Transfect to (mL)` = volume_transfect,
                `Transfect\nto` = plate_count
            ) %>%
            mutate(across(contains("(µL)"), round_volumes)) %>%
            select(-`Transfect to (mL)`)
        
        # Function to calculate master mix
        calc_master_mix <- function(sample_data, plate_input) {
            
            # safely calculate total for each reagent column
            calc_total <- function(col) {
                sample_data() %>%
                    pull(col) %>%
                    sum(na.rm = TRUE)
            }
            
            # totals
            total_tf_vol <- calc_total("volume_transfect")
            total_hibit_vol <- calc_total("volume_hibit")
            total_luc2_vol <- calc_total("volume_luc2")
            
            # total plates (avoid divide-by-zero)
            tot_plates <- ifelse(
                as.numeric(plate_input) == 0, 
                0,
                total_tf_vol / as.numeric(plate_input)
            )
            
            # helper for master volume
            calc_master <- function(total_vol) {
                mm_extra <- ifelse(input$mm_extra < 0, 0, input$mm_extra)
                
                ifelse(
                    tot_plates == 0, 
                    0, 
                    (total_vol / tot_plates) * (tot_plates + mm_extra)
                )
            }
            
            # return named list of results
            list(
                total_tf_vol = total_tf_vol,
                total_hibit_vol = total_hibit_vol,
                total_luc2_vol = total_luc2_vol,
                tot_plates = tot_plates,
                hibit_master_vol = calc_master(total_hibit_vol),
                luc2_master_vol = calc_master(total_luc2_vol)
            )
        }
        
        # Generate master mix text only if needed
        if (input$show_mm) {
            results <- calc_master_mix(sample_data, input$plate_input)
            mastermix_text <- paste(
                "To create a HiBiT + Luc2 master mix, add",
                round_volumes(results$hibit_master_vol),
                "uL HiBiT plasmid and",
                round_volumes(results$luc2_master_vol),
                "uL Luc2 plasmid to a master mix tube,",
                "then add the indicated master mix volume to each tube.\n"
            )
        } else {
            data <- select(data, -`Master mix\n(µL)`)
        }
        
        # Remove HiBiT and Luc2 volumes if turned off in settings
        if (!input$show_hibit_luc2) {
            data <- select(data, -`HiBiT\n(µL)`, -`Luc2\n(µL)`)
        }
        
        # Protocol text
        protocol_text <- paste(
            "Protocol:\n",
            "1. Add calculated volumes of DNA to a sterile tube.\n",
            "2. Add calculated volume of Opti-MEM.\n",
            "3. Add calculated volume of Trans-IT.\n",
            "4. Incubate samples for 15 minutes at room temperature.\n",
            "5. Add sample volume equivalent to 10% of cell medium volume.", 
            "(For example, add 2 mL to each 20 mL dish.)"
        )
        
        # process into flextable
        set_flextable_defaults(
            font.family = "Helvetica",
            font.size = 14,
            word_wrap = FALSE
        )
        
        flextable <- data %>%
            flextable() %>%
            autofit() %>%
            add_header_lines(
                values = paste0(
                    "Pseudovirus transfection table (", Sys.Date(), ")"
                )
            ) %>%
            add_footer_lines(
                values = as_paragraph(
                    paste("HiBiT concentration:", input$hibit_input, "ng/µL"),
                    "\n",
                    paste("Luc2 concentration:", input$luc2_input, "ng/µL")
                )
            )
        
        # conditionally add mastermix_text
        if (input$show_mm) {
            flextable <- flextable %>%
                add_footer_lines(values = as_paragraph(mastermix_text))
        }
        
        # always add protocol_text
        flextable <- flextable %>%
            add_footer_lines(values = as_paragraph(protocol_text)) %>%
            bold(bold = TRUE, part = "header") %>%
            bold(j = 3:6, bold = TRUE, part = "body") %>%
            color(j = 3:5, color = "red", part = "body") %>%
            color(j = 6, color = "blue", part = "body") %>%
            align(align = "center", part = "header") %>%
            align(align = "center", part = "body")
        
        flextable
    })
    
    # Button - add sample
    observeEvent(input$add_sample, {
        # calculate vols and consolidate into new row
        vols <- calc_volumes()
        new_row <- data.frame(
            sample_id = input$sample_name,
            conc_spike = input$spike_input,
            volume_spike = vols$volume_spike_calc,
            volume_hibit = vols$volume_hibit_calc,
            volume_luc2 = vols$volume_luc2_calc,
            volume_master = vols$volume_master_calc,
            volume_optimem = vols$volume_optimem_calc,
            volume_transit = vols$volume_transit_calc,
            volume_transfect = vols$volume_cell_medium,
            plate_count = vols$plate_count
        )
        
        # bind to current
        current <- sample_data()
        sample_data(rbind(current, new_row))
    })
    
    # button - remove all samples by setting sample_data to empty.
    observeEvent(input$remove_all_samples, {
        sample_data(empty_table)
    })
    
    # Reactive re-calculation of hibit & luc2 volumes when inputs change
    observeEvent({
        list(input$hibit_input, input$luc2_input)
    }, {
        # Only recalc if there’s existing data
        if (nrow(sample_data()) > 0) {
            df <- sample_data()
            vols <- calc_volumes()
            
            # Update hibit, luc2, and master mix volumes in existing table
            df$volume_hibit <- (400 * df$volume_transfect) / input$hibit_input
            df$volume_luc2  <- (400 * df$volume_transfect) / input$luc2_input
            df$volume_master <- df$volume_hibit + df$volume_luc2
            
            sample_data(df)
        }
    })
    
    # render flextable as HTML widget
    output$table_output <- renderUI({
        validate_inputs()
        
        tags$div(
            style = "margin-left: 0; margin-right: auto; width: fit-content;",
            flextable::htmltools_value(ft())
        )
    })

    # render flextable as .docx and make downloadable
    output$download_docx <- downloadHandler(
        filename = function() {
            paste0("transfection_table_", Sys.Date(), ".docx")
        },
        content = function(file) {
            # format for docx. Increment integer to increase width
            ft_docx <- ft() %>%
                width(
                    width = dim(.)$widths * 10.5 / (flextable_dim(.)$widths)
                )

            # format to A4 landscape with narrow margins
            sect_properties <- prop_section(
                page_size = page_size(
                    orient = "landscape",
                    width = 11.69, 
                    height = 8.27
                ),
                type = "continuous",
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
