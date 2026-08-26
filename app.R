suppressPackageStartupMessages({
  library(pedtools)
  library(pedprobr)
  library(pedmut)
  library(pedFamilias)
  library(shiny)
  library(bslib)
})

ui = page_sidebar(
  tags$head(
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$script(src = "scripts.js")
  ),
  
  title = div(
    class = "app-title",
    tags$img(src = "pedlogo.svg", class = "app-logo"),
    div(
      class = "app-title-text",
      div(class = "app-name", "famnesia"),
      div(class = "app-subtitle", "Anonymising Familias files")
    )
  ),
  #theme = bs_theme(version = 5, primary = "#3f7f80", navbar_bg = "#eee8df"),
  theme = bs_theme(version = 5, primary = "#526f8e", navbar_bg = "#e7e5e1"),
  fillable = FALSE,

  sidebar = sidebar(
    width = 260,
    open = "always",
    resizable = FALSE,
    gap = "1.5rem",
    bg = "#f8f7f4",


    fileInput("fileInput", "Familias .fam file", buttonLabel = icon("folder-open"), accept = ".fam"),
    div(
      class = "file-actions",
      actionButton("example", "Example",
                   class = "btn-sm btn-outline-secondary"),
      actionButton("settings", "Settings", icon = icon("gear"),
                   class = "btn-sm btn-outline-secondary")
    ),
    checkboxGroupInput(
      "options", "Masking options",
      choiceNames = list(
        "Family names"  |> addTip("Rename families to F1,F2,..."),
        "ID labels"     |> addTip("Rename individuals to 1,2,..."),
        "Marker names"  |> addTip("Rename markers randomly to M1, M2, ..."),
        "Lump alleles"  |> addTip("Merge unobserved alleles at each marker"),
        "Randomize sex" |> addTip("Swap sex of random (suitable) individuals")
      ),
      choiceValues = c("famnames", "ids", "markernames", "lump", "sex"),
      selected = c("famnames", "ids", "markernames")
    ),
    radioButtons("alleles", "Allele labels",
      choiceNames = list(
        "1,2,3,..." |> addTip("Rename alleles randomly to 1,2,..."),
        "Constrained" |> addTip("Rename, but preserve order and decimal groups")
      ),
      choiceValues = c("strong", "constrained")
    ),
    radioButtons("freqs", "Frequencies",
      choiceNames = list(
        "Original" |> addTip("Keep original freqs"),
        "Tweak"    |> addTip("Perturb up to 2% and renormalize"),
        "Round"    |> addTip("Round to 4 decimals")
      ),
      choiceValues = c("original", "tweak", "round")
    ),
    radioButtons("mutmodels", "Mutation models",
      choiceNames = list(
        "Original" |> addTip("Keep original parameters"),
        "Simplify" |> addTip("Use 'equal' model with original rates"),
        "Disable"  |> addTip("Remove all models")
      ),
      choiceValues = c("original", "simplify", "disable")
    ),
    div(
      class = "d-flex gap-2",
      actionButton("apply", "Apply", class = "btn-sm btn-primary flex-fill"),
      downloadButton("download", "Download",
                     class = "btn-sm btn-outline-secondary flex-fill text-nowrap")
    )
  ),

  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header(
        div(class = "card-heading",
            span("Original  ", textOutput("statsOriginal", inline = TRUE)),
            uiOutput("lrOriginal", inline = TRUE)),
        class = "original-heading"
      ),
      DT::DTOutput("tableOriginal", fill = FALSE),
      plotOutput("plotOriginal", height = "320px")
    ),
    card(
      card_header(
        div(class = "card-heading",
            span("Masked  ", textOutput("statsMasked", inline = TRUE)),
            uiOutput("lrMasked", inline = TRUE)),
        class = "masked-heading"
      ),
      DT::DTOutput("tableMasked", fill = FALSE),
      plotOutput("plotMasked", height = "320px")
    )
  )
)

server = function(input, output, session) {
  imported = reactiveVal(NULL)
  masked = reactiveVal(NULL)
  prefs = reactiveValues(removeEmpty = TRUE,
                         abbreviate = TRUE,
                         noMutLR = FALSE)


  # Import data ---------------------------------------------------------------------------------
  
  # Shared pipeline for upload & example
  importFam = function(path) {
    imported(NULL)
    masked(NULL)
  
    tryCatch({
      x = loadFamData(path)
      imported(x)
  
      if((x$params$theta %||% 0) > 0)
        showNotification("Note: theta ignored in LR calculations", duration = 5)
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  }
  
  observeEvent(input$fileInput, importFam(req(input$fileInput$datapath)))
  observeEvent(input$example, importFam("data/sibship.fam"))


  # Apply masking -------------------------------------------------------------------------------
  
  # Working version of input: Remove empty markers if indicated
  original = reactive({
    x = req(imported())
    keep = !prefs$removeEmpty | lengths(x$observed) > 0
    markers = names(x$attrs)[keep]

    if(!all(keep))
      x$peds = lapply(x$peds, selectMarkers, markers = markers)

    x$attrs = x$attrs[keep]
    x$observed = x$observed[keep]
    x$mutpars = x$mutpars[keep]
    x$lr = x$lr[keep]
    if(!is.null(x$lrNoMut))
      x$lrNoMut = x$lrNoMut[keep]
    x
  })

  # Build a fresh masked copy from the converted data
  observeEvent(input$apply, {
    masked(NULL)
    options = input$options
    
    tryCatch({
      mdat = maskData(
        req(original()),
        famnames = "famnames" %in% options,
        ids = "ids" %in% options,
        markernames = "markernames" %in% options,
        lump = "lump" %in% options,
        sex = "sex" %in% options,
        alleles = input$alleles,
        freqs = input$freqs,
        mutmodels = input$mutmodels
      )
      masked(mdat)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })


  # Main tables and plots -----------------------------------------------------------------------
  
  output$tableOriginal = DT::renderDT({
    input$apply
    x = req(original())
    renderMarkerTable(x$peds, x$attrs, x$lr, 
                      lrNoMut = if(prefs$noMutLR) x$lrNoMut,
                      shortNames = prefs$abbreviate)
  }, server = FALSE)

  output$tableMasked = DT::renderDT({
    x = req(original())
    m = req(masked())
    prefs$abbreviate # Refresh both tables together
    renderMarkerTable(m$peds, m$attrs, m$lr,
                      lrNoMut = if(prefs$noMutLR) m$lrNoMut,
                      referenceLR = unname(x$lr))
  }, server = FALSE)

  output$statsOriginal = renderText(pedStats(req(original())))
  output$statsMasked = renderText(pedStats(req(masked())))
  
  output$lrOriginal = renderUI({
    lr = prod(req(original())$lr)
    tags$span(class = "lr-total",  paste0("LR=", sprintf("%.3g", lr)))
  })

  output$lrMasked = renderUI({
    lr0 = prod(req(original())$lr)
    lr = prod(req(masked())$lr)
    dev = 100 * (lr / lr0 - 1)
    cls = if(abs(dev) < 1) "lr-close" else if(abs(dev) < 5) "lr-medium" else "lr-large"

    tags$span(
      class = "lr-summary",
      tags$span(class = "lr-total", paste0("LR=", sprintf("%.3g", lr))),
      tags$span(class = paste("lr-change", cls), sprintf("%+.1f%%", dev))
    )
  })

  # Pedigree plots
  output$plotOriginal = renderPlot({
    plotAllPeds(req(original())$peds)
  }, execOnResize = TRUE)

  output$plotMasked = renderPlot({
    plotAllPeds(req(masked())$peds)
  }, execOnResize = TRUE)


  # Frequency modal -----------------------------------------------------------------------------
  
  markerDetails = reactive({
    i = req(input$markerDblclick)
    x = req(original())
    m = req(masked())

    list(original = x$attrs[[i]],
         masked = m$attrs[[i]],
         alleleMap = m$alleleMap[[i]],
         markerNames = c(names(x$attrs)[i], names(m$attrs)[i]))
  })
  
  observeEvent(input$markerDblclick, {
    d = req(markerDetails())

    showModal(modalDialog(
      title = tags$span(
        class = "freq-title",
        tags$span(paste("Alleles and frequencies:",
                        paste(d$markerNames, collapse = " → "))),
        modalButton("Close")
      ),
      DT::DTOutput("freqTable"),
      size = "m",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Frequency table (in modal)
  output$freqTable = DT::renderDT({
    d = markerDetails()
    renderFreqTable(d$original, d$masked, d$alleleMap)
  }, server = FALSE)
      

  # Download ------------------------------------------------------------------------------------
        
  output$download = downloadHandler(
    filename = function() "famnesia.fam",
    content = function(file) {
      m = req(masked())
      famfile = tempfile(fileext = ".fam")
      on.exit(unlink(famfile), add = TRUE)

      writeFam(m$peds, famfile = famfile, params = m$params, verbose = FALSE)

      if(!file.copy(famfile, file, overwrite = TRUE))
        stop("Could not prepare the download", call. = FALSE)
    }
  )


  # Settings ------------------------------------------------------------------------------------

  observeEvent(input$settings, {
    showModal(modalDialog(
      title = tagList(icon("gear"), "Settings"),
      tags$div(
        class = "border rounded-3 bg-body-tertiary p-3",
        tags$h6("Data and display", class = "mb-3"),
        checkboxInput(
          "settingRemoveEmpty", "Remove empty markers",
          value = prefs$removeEmpty
        ),
        checkboxInput(
          "settingAbbreviate", "Abbreviate long names in Original table",
          value = prefs$abbreviate
        ),
        checkboxInput(
          "settingNoMutLR", "Include LRs without mutation models",
          value = prefs$noMutLR
        )
      ),
      size = "s",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$settingRemoveEmpty, {
    value = input$settingRemoveEmpty
    if(!identical(value, prefs$removeEmpty)) {
      prefs$removeEmpty = value
      masked(NULL)
    }
  })
  
  observeEvent(input$settingAbbreviate, {
    if(!identical(input$settingAbbreviate, prefs$abbreviate))
      prefs$abbreviate = input$settingAbbreviate
  })
  
  observeEvent(input$settingNoMutLR, {
    if(!identical(input$settingNoMutLR, prefs$noMutLR))
      prefs$noMutLR = input$settingNoMutLR
  })

}

shinyApp(ui = ui, server = server)
