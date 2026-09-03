suppressPackageStartupMessages({
  library(pedtools)
  library(pedprobr)
  library(pedmut)
  library(pedFamilias)
  library(shiny)
  library(bslib)
})

ui = page_sidebar(
  includeCSS("www/styles.css"),
  tags$head(tags$script(src = "scripts.js")),
  shinyjs::useShinyjs(),
  useBusyIndicators(),
  
  title = div(
    class = "app-title",
    tags$img(src = "pedlogo.svg", class = "app-logo"),
    div(
      class = "app-title-text",
      div(class = "app-name", "Famnesia"),
      div(class = "app-subtitle", "Anonymising Familias files")
    ),
    actionLink("settings",
           tagList(icon("gear", class = "fs-3"), tags$br(), tags$small("Settings")),
           class = "app-settings text-center text-secondary")
  ),
  theme = bs_theme(version = 5, primary = "#526f8e", 
                   navbar_bg = "#e7e5e1"),
  fillable = FALSE,

  sidebar = sidebar(
    width = 260,
    open = "always",
    resizable = FALSE,
    gap = "0.9rem",
    bg = "#f8f7f4",

    div(
      class = "file-input-wrap",
      div(
        class = "d-flex align-items-center justify-content-between fw-semibold",
        span("Familias file"),
        actionButton("example", "EXAMPLE", class = "btn-sm example-btn")
      ),
      fileInput("fileInput", NULL, buttonLabel = icon("folder-open"),
                accept = ".fam")
    ),
    actionButton("analyse", "Analyse", icon = icon("magnifying-glass"),
             class = "btn-sm btn-outline-primary"),
    
    checkboxGroupInput(
      "options", "Masking options",
      choiceNames = list(
        "Family names"  |> addTip("Rename families to F1,F2,..."),
        "ID labels"     |> addTip("Rename individuals to 1,2,..."),
        "Marker names"  |> addTip("Rename markers to M1,M2,..."),
        "Shuffle markers"  |> addTip("Permute the marker order"),
        "Lump alleles"  |> addTip("Merge unobserved alleles at each marker"),
        "Randomize sex" |> addTip("Swap sex of random (suitable) individuals")
      ),
      choiceValues = c("famnames", "ids", "markernames", "shuffle", "lump", "sex"),
      selected = c("famnames", "ids", "markernames", "shuffle")
    ),
    radioButtons("alleles", "Allele labels",
      choiceNames = list(
        "Original" |> addTip("Keep original allele labels"),
        "Constrained" |> addTip("Rename, but preserve order and decimal groups"),
        "1,2,3,..." |> addTip("Rename alleles randomly to 1,2,...")
      ),
      choiceValues = c("original", "constrained", "strong")
    ),
    radioButtons("freqs", "Frequencies",
      choiceNames = list(
        "Original" |> addTip("Keep original freqs"),
        "Round"    |> addTip("Round to 4 decimals"),
        "Tweak"    |> addTip("Perturb up to 2% and renormalize")
      ),
      choiceValues = c("original", "round", "tweak")
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
      class = "d-flex gap-1",
      actionButton("apply", "Apply", icon = icon("wand-magic-sparkles"),
                   class = "btn-primary flex-fill text-nowrap px-2"),
      downloadButton("download", "Download",
                     class = "btn-outline-secondary flex-fill text-nowrap px-2")
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
      plotOutput("plotOriginal", fill = FALSE, height = "320px")
    ),
    card(
      card_header(
        div(class = "card-heading",
            span("Masked  ", textOutput("statsMasked", inline = TRUE)),
            uiOutput("lrMasked", inline = TRUE)),
        class = "masked-heading"
      ),
      div(style = "position: relative;",
        DT::DTOutput("tableMasked", fill = FALSE),
        uiOutput("maskedHint")
      ),
      plotOutput("plotMasked", fill = FALSE, height = "320px")
    )
  )
)

server = function(input, output, session) {
  
  if(interactive())
    session$onSessionEnded(stopApp)
  
  # Main reactives
  imported = reactiveVal(NULL)
  masked = reactiveVal(NULL)
  prefs = reactiveValues(removeEmptyComps = TRUE,
                         removeEmptyMarkers = TRUE,
                         abbreviate = TRUE,
                         noMutLR = FALSE,
                         lrDigits = 2,
                         seed = 12345)


  # Import data ---------------------------------------------------------------------------------
  
  # Shared pipeline for upload & example
  importFam = function(path, filename = basename(path)) {
    imported(NULL)
    masked(NULL)
  
    tryCatch({
      x = loadFamData(path)
      x$filename = filename
      imported(x)
  
      if((x$params$theta %||% 0) > 0)
        showNotification("Note: theta ignored in LR calculations", duration = 5)
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 5)
    })
  }
  
  observeEvent(input$fileInput, {d = req(input$fileInput); importFam(d$datapath, d$name)})
  
  observeEvent(input$example, importFam("data/sibship.fam"))

  
  setMasking = function(choices) {
    nms = names(choices)
    if("options" %in% nms)
      updateCheckboxGroupInput(session, "options", selected = choices$options)
    if("alleles" %in% nms)
      updateRadioButtons(session, "alleles", selected = choices$alleles)
    if("freqs" %in% nms)
      updateRadioButtons(session, "freqs", selected = choices$freqs)
    if("mutmodels" %in% nms)
      updateRadioButtons(session, "mutmodels", selected = choices$mutmodels)
  }

  # File analysis -----------------------------------------------------------------------------
  
  analysis = reactive(analyseFam(req(original())))
  
  observeEvent(input$analyse, {
    a = analysis()
  
    row = function(id, label, text, effect, rec = FALSE, cls = "btn-outline-dark")
      div(class = paste("row align-items-center px-2 py-2 border-top text-dark",
                  if(rec) "bg-success bg-opacity-10"),
          div(class = "col-3",
              actionButton(id, label, width = "100%",
                           class = paste(cls, "text-center text-nowrap px-2"))),
          div(class = "col-7 small", text),
          div(class = "col-2 text-center fs-3", effectSymbol(effect)))
  
    showModal(modalDialog(
      title = tags$span(
        class = "freq-title",
        tagList(icon("magnifying-glass"), "File analysis"),
        modalButton("Close")
      ),
      
      h5("File summary"),
      
      div(class = "bg-light rounded-3 p-3 mb-3",
          tags$p(class = "mb-0",
            tags$b("File name: "), a$info$filename, tags$br(),
            tags$b("Family names: "), paste(a$info$famnames, collapse = ", "), tags$br(),
            tags$b("Typed indivs: "), a$info$individuals, tags$br(),
            tags$b("Database: "), basename(a$info$database), tags$br(),
            tags$b("Mutation models: "), a$info$mutmodel
      )),
      
      tags$hr(class = "my-3"),
      h5("Masking strategies for this file"),
      
      div(class = "border rounded-3 overflow-hidden mb-2",
          tags$div(class = "row px-2 py-2 fw-semibold bg-light",
            div(class = "col-3"),
            div(class = "col-7", "Settings"),
            div(class = "col-2 text-center", "LR")
          ),
          row("analyseMax", "Max", maskingText(a$max), "non-exact"),
          row("analyseSuggested", "Suggested", maskingText(a$suggested), a$effect, rec = TRUE, 
              cls = "btn-success"),
          row("analysePreserve", "Keep LR", maskingText(a$preserve), "exact", 
              cls = "btn-outline-primary")
      ),
      tags$small(class = "text-secondary", "Standard masks are included in all cases."),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
    
  observeEvent(input$analyseMax, {setMasking(analysis()$max); removeModal()})
  observeEvent(input$analyseSuggested, {setMasking(analysis()$suggested); removeModal()})
  observeEvent(input$analysePreserve, {setMasking(analysis()$preserve); removeModal()})
  
  
  observe({
    orig = req(original())
    stepwise = any(vapply(orig$mutpars, FUN.VALUE = logical(1),
                          \(p) !is.null(p) && "stepwise" %in% unlist(p$model)))
    incompatible = stepwise && "lump" %in% input$options
  
    shinyjs::toggleState(
      selector = "#mutmodels input[value='original']",
      condition = !incompatible)
  
    if(incompatible && input$mutmodels == "original")
      updateRadioButtons(session, "mutmodels", selected = "simplify")
  })
  
  # Apply masking -------------------------------------------------------------------------------
  
  # Working version of input: Remove empty markers if indicated
  original = reactive({
    x = req(imported())
    keep = !prefs$removeEmptyMarkers | lengths(x$observed) > 0
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
    orig = req(original())
    masked(NULL)
    options = input$options
    if(!is.null(prefs$seed)) 
      set.seed(prefs$seed)
    
    tryCatch({
      mdat = maskData(
        orig,
        famnames = "famnames" %in% options,
        ids = "ids" %in% options,
        markernames = "markernames" %in% options,
        shuffle = "shuffle" %in% options,
        lump = "lump" %in% options,
        sex = "sex" %in% options,
        alleles = input$alleles,
        freqs = input$freqs,
        mutmodels = input$mutmodels
      )
      masked(mdat)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 5)
    })
  })


  # Main tables and plots -----------------------------------------------------------------------
  
  output$tableOriginal = DT::renderDT({
    input$apply
    x = req(original())
    renderMarkerTable(x$peds, x$attrs, x$lr, digits = prefs$lrDigits,
                      shortNames = prefs$abbreviate,
                      lrNoMut = if(prefs$noMutLR) x$lrNoMut)
  }, server = FALSE)

  output$tableMasked = DT::renderDT({
    x = req(original())
    m = req(masked())
    renderMarkerTable(m$peds, m$attrs, m$lr, digits = prefs$lrDigits,
                      shortNames = prefs$abbreviate,
                      lrNoMut = if(prefs$noMutLR) m$lrNoMut,
                      referenceLR = unname(x$lr))
  }, server = FALSE)

  output$statsOriginal = renderText(pedStats(req(original())))
  output$statsMasked = renderText(pedStats(req(masked())))
  
  output$lrOriginal = renderUI({
    orig = req(original())
    LRtag(prod(orig$lr), length(orig$peds), prefs$lrDigits)
  })
  
  output$lrMasked = renderUI({
    orig = req(original())
    npeds = length(orig$peds)
    lr0 = prod(orig$lr)
    lr = prod(req(masked())$lr)
  
    if(!is.finite(lr0) || !is.finite(lr) || lr0 == 0)
      return(LRtag(lr, npeds, prefs$lrDigits))
    
    dev = 100 * (lr / lr0 - 1)
    cls = if(abs(dev) < 1) "lr-close" else if(abs(dev) < 5) "lr-medium" else "lr-large"
    
    tags$span(
      class = "lr-summary",
      LRtag(lr, npeds, prefs$lrDigits),
      tags$span(class = paste("lr-change", cls), sprintf("%+.1f%%", dev))
    )
  })
  
  output$maskedHint = renderUI({
    if(is.null(masked())) return(NULL)
    div(class = "table-hint", "Double-click a row to view alleles and frequencies.")
  })
    
  # Pedigree plots
  output$plotOriginal = renderPlot({
    plotAllPeds(req(original())$peds, removeEmpty = prefs$removeEmptyComps)
  }, execOnResize = TRUE)

  output$plotMasked = renderPlot({
    plotAllPeds(req(masked())$peds, removeEmpty = prefs$removeEmptyComps)
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

      peds = lapply(m$peds, selectMarkers, markers = m$markerOrder)
      writeFam(peds, famfile = famfile, params = m$params, verbose = FALSE)
      
      if(!file.copy(famfile, file, overwrite = TRUE))
        stop("Could not prepare the download", call. = FALSE)
    }
  )


  # Settings ------------------------------------------------------------------------------------

  observeEvent(input$settings, {
    showModal(modalDialog(
      title = tagList(icon("gear"), "Settings"),
      tags$div(
        class = "border rounded-3 bg-body-tertiary p-3 text-nowrap",
        checkboxInput("settingRemoveEmptyMarkers", "Remove empty markers",
                      value = prefs$removeEmptyMarkers),
        checkboxInput("settingRemoveEmptyComps", "Hide empty plot components",
                      value = prefs$removeEmptyComps),
        checkboxInput("settingAbbreviate", "Shorten long names in table",
                      value = prefs$abbreviate),
        checkboxInput("settingNoMutLR", "Include LRs without mutation models",
                      value = prefs$noMutLR),
        numericInput("settingLRdigits", "LR decimals", min = 0, max = 10, step = 1,
                     value = prefs$lrDigits),
        numericInput("settingSeed", "Random seed", min = 0, step = 1, 
                     value = prefs$seed)
      ),
      size = "m",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$settingRemoveEmptyComps, {
    value = input$settingRemoveEmptyComps
    if(!identical(value, prefs$removeEmptyComps)) {
      prefs$removeEmptyComps = value
    }
  })
  
  observeEvent(input$settingRemoveEmptyMarkers, {
    value = input$settingRemoveEmptyMarkers
    if(!identical(value, prefs$removeEmptyMarkers)) {
      prefs$removeEmptyMarkers = value
      masked(NULL)
    }
  })
  
  observeEvent(input$settingAbbreviate, {
    if(!identical(input$settingAbbreviate, prefs$abbreviate))
      prefs$abbreviate = input$settingAbbreviate
  })
  
  observeEvent(input$settingLRdigits, {
    d = input$settingLRdigits
    if(is.finite(d)) 
      prefs$lrDigits = max(0, min(10, round(d)))
  })
  
  observeEvent(input$settingSeed, {
    s = input$settingSeed
    prefs$seed = if(is.na(s)) NULL else s
  })
  
  observeEvent(input$settingNoMutLR, {
    if(!identical(input$settingNoMutLR, prefs$noMutLR))
      prefs$noMutLR = input$settingNoMutLR
  })

}

shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
