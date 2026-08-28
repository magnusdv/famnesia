`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

setnames = function(x, nms = x) {
  names(x) = nms
  x
}

uniqueAlleles = function(x)
  x[!is.na(x)] |> as.character() |> unique.default()

formatLR = function(z, digits = 2)
  ifelse(is.na(z), NA_character_, sprintf(paste0("%.", digits, "f"), z))

LRtag = function(lr, npeds, digits) {
    fmt = paste0("%.", digits + 1, "g")
    tags$span("LR", if(npeds > 2) tags$sub("1:2"), " = ", sprintf(fmt, lr),
              class = "lr-total")
}

# Truncated allele list for table display
shortAlleles = function(x, show = 5, maxlen = 20) {
  x = as.character(x)
  if(length(x) > show)
    x = c(x[1:(show - 2)], "...", tail(x, 1))
  s = paste(x, collapse = ", ")
  if(nchar(s) <= maxlen || show <= 3)
    s
  else
    shortAlleles(x, show - 1, maxlen)
}

pedStats = function(x)
  sprintf("(%d peds, %d typed, %d markers)",
          length(x$peds), length(typedMembers(x$peds[[1]])), length(x$attrs))

addTip = function(label, text)
  tooltip(tags$span(label), text, placement = "right")



loadFamData = function(path) {
  famfile = tempfile(fileext = ".fam")
  on.exit(unlink(famfile), add = TRUE)

  copied = file.copy(path, famfile)
  if(!copied)
    stop("Could not read FAM file", call. = FALSE)

  x = readFam(famfile, convert = FALSE, includeParams = TRUE, verbose = FALSE)

  if(isTRUE(x$params$dvi))
    stop("DVI files are not supported", call. = FALSE)
  if(!length(x$pedigrees))
    stop("The file contains no pedigrees", call. = FALSE)

  orig = Familias2ped(x$pedigrees, datamatrix = x$datamatrix, loci = x$loci)
  
  if(inherits(x$pedigrees, "FamiliasPedigree"))
    orig = list(orig)

  attrs = getLocusAttributes(orig[[1]], attribs = c("alleles", "afreq", "mutmod"))
  markers = names(attrs)
  Nm = length(markers)
  
  # Read observed alleles from datamatrix
  if(is.null(x$datamatrix))
    observed = rep(list(character()), Nm)
  else
    observed = lapply(seq_len(Nm), \(i) uniqueAlleles(x$datamatrix[, 2*i - 1:0]))
  
  names(observed) = markers

  # Save mutmodel parameters
  mutpars = lapply(attrs, \(a)
    if(is.null(a$mutmod)) NULL else getParams(a$mutmod, format = 4))
  
  lr = markerLR(orig, theta = 0)
  lrNoMut = if(any(lengths(mutpars)))
    markerLR(lapply(orig, setMutmod, model = NULL), theta = 0)

  list(
    peds = orig,
    attrs = attrs,
    observed = observed,
    mutpars = mutpars,
    lr = lr,
    lrNoMut = lrNoMut,
    params = x$params
  )
}
   
plotAllPeds = function(peds, removeEmpty = TRUE) {
  if(removeEmpty)
    peds = removeEmptyComps(peds)
  npeds = length(peds)
  cex = if(npeds == 1) 1.1 else 1.3

  tryCatch(
    suppressWarnings(plotPedList(peds, hatched = typedMembers, cex = cex)),
    error = function(e) {
      msg = conditionMessage(e)
      if(grepl("Cannot fit the graph", msg))
        msg = "Pedigrees do not fit.\nTry enlarging the window."
      shiny::validate(msg, errorClass = "plot")
    }
  )
}

removeEmptyComps = function(x) {
  if(is.null(x)) 
    return(NULL)
  for(pedname in names(x)) {
    ped = x[[pedname]]
    if(is.ped(ped))
      next
    empty = vapply(ped, function(comp) !any(unlist(comp$MARKERS)), FALSE)
    if(all(empty)) 
      warning(sprintf("Pedigree '%s' has no typed members", pedname))
    else x[[pedname]] = ped[!empty]
  }
  x
}
    

# Frequency table (in modal)

renderFreqTable = function(orig, mask, map) {
  new = unname(map[orig$alleles])
  tab = data.frame(
    orig$alleles,
    unname(orig$afreq),
    new,
    unname(mask$afreq[match(new, mask$alleles)])
  )
  tab[[3]][is.na(tab[[3]])] = ""

  # Add the lump as a separate row
  lump = unname(map["lump"])
  if(!is.na(lump))
    tab = rbind(tab, list("", NA_real_, paste(lump, "(lump)"),
                          unname(mask$afreq[match(lump, mask$alleles)])))

  tab[c(2, 4)] = lapply(tab[c(2, 4)], signif, 6)

  header = htmltools::HTML(
    "<table><thead>
     <tr>
      <th colspan='2' class='freq-group'>Original</th>
      <th colspan='2' class='freq-group masked'>Masked</th>
     </tr>
     <tr>
      <th>Allele</th><th>Frequency</th>
      <th class='masked'>Allele</th><th>Frequency</th>
     </tr>
     </thead></table>"
  )

  DT::datatable(
    tab, container = header, rownames = FALSE,
    class = "compact stripe",
    options = list(
      dom = "t", paging = FALSE, ordering = TRUE,
      scrollY = "55vh", scrollCollapse = TRUE,
      orderCellsTop = FALSE,
      columnDefs = list(list(
        targets = 0:3, className = "dt-center",
        render = DT::JS("sortBlanksLast")
      ))
    )
  ) |> 
    DT::formatStyle(1:4,  lineHeight = "98%")
}
