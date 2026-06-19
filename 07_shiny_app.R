################################################################################
# app.R
################################################################################

library(shiny)
library(glmnet)
library(dplyr)

clock_models <- readRDS("results/models/clock_models.rds")
clock_calibration <- readRDS("results/models/clock_calibration.rds")

model_choices <- c(
  "Inflammation-MR-Clock1" = "Inflammation_MR_Clock1",
  "Inflammation-MR-Clock2" = "Inflammation_MR_Clock2",
  "Proteomics-MR-Clock1" = "Proteomics_MR_Clock1",
  "Proteomics-MR-Clock2" = "Proteomics_MR_Clock2"
)

get_model_key <- function(clock, sex){
  paste(clock, sex, sep = "_")
}

get_model_variables <- function(clock, sex){
  key <- get_model_key(clock, sex)
  clock_models[[key]]$variables
}

predict_mr_age <- function(input_values, clock, sex){
  key <- get_model_key(clock, sex)
  model_obj <- clock_models[[key]]
  calibration_obj <- clock_calibration[[key]]

  vars <- model_obj$variables
  x <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(vars)))
  names(x) <- vars

  for(v in vars){
    x[[v]] <- as.numeric(input_values[[v]])
  }

  x_mat <- as.matrix(x)

  raw_age <- as.numeric(
    predict(model_obj$fit,
            newx = x_mat,
            s = "lambda.min")
  )

  mr_age <- as.numeric(
    predict(calibration_obj,
            newdata = data.frame(pred_train_raw = raw_age))
  )

  mr_age
}

ui <- fluidPage(
  titlePanel("Inflammation- and Proteomics-MR Ageing Clock Calculator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("sex",
                  "Sex",
                  choices = c("Men", "Women")),
      selectInput("clock",
                  "MR ageing clock",
                  choices = model_choices),
      numericInput("age",
                   "Chronological age",
                   value = 60,
                   min = 18,
                   max = 100),
      uiOutput("dynamic_inputs"),
      actionButton("calculate",
                   "Calculate MR age")
    ),
    mainPanel(
      h3("Estimated biological age"),
      verbatimTextOutput("clock_output"),
      helpText("MR age is estimated using sex-specific elastic-net ageing clocks and linear calibration models derived in the UK Biobank derivation set.")
    )
  )
)

server <- function(input, output, session){

  output$dynamic_inputs <- renderUI({
    req(input$clock, input$sex)

    vars <- get_model_variables(input$clock, input$sex)

    lapply(vars, function(v){
      numericInput(inputId = v,
                   label = v,
                   value = 0,
                   step = 0.01)
    })
  })

  result <- eventReactive(input$calculate, {
    req(input$clock, input$sex, input$age)

    vars <- get_model_variables(input$clock, input$sex)

    input_values <- lapply(vars, function(v){
      input[[v]]
    })
    names(input_values) <- vars

    mr_age <- predict_mr_age(input_values,
                             input$clock,
                             input$sex)

    acceleration <- mr_age - input$age

    list(
      mr_age = mr_age,
      acceleration = acceleration
    )
  })

  output$clock_output <- renderPrint({
    req(result())

    cat("Selected clock:", names(model_choices)[model_choices == input$clock], "\n")
    cat("Sex:", input$sex, "\n")
    cat("Chronological age:", round(input$age, 1), "years\n")
    cat("Estimated MR age:", round(result()$mr_age, 1), "years\n")
    cat("Age acceleration:", round(result()$acceleration, 1), "years\n")
  })
}

shinyApp(ui, server)
