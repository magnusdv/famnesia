`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

setnames = function(x, nms = x) {
  names(x) = nms
  x
}

uniqueAlleles = function(x)
  x[!is.na(x)] |> as.character() |> unique.default()

formatLR = function(z)
  ifelse(is.na(z), NA_character_, sprintf("%.2f", z))

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
   
      
# Marker-wise LRs
markerLR = function(peds, theta = 0) {
  lik = do.call(cbind, lapply(peds, pedprobr::likelihood, theta = theta))
  lr = if(ncol(lik) < 2) rep(NA_real_, nrow(lik)) else lik[, 1] / lik[, 2]
  setnames(lr, name(peds[[1]]))
}

renderMarkerTable = function(peds, locusAttributes, lr, lrNoMut = NULL,
                             referenceLR = NULL, shortNames = FALSE) {
  markers = names(locusAttributes)
  attrs = locusAttributes[markers]
  
  # Prepare genotype columns
  p1 = peds[[1]]
  geno = getGenotypes(p1, typedMembers(p1), markers = markers) |> t.default()
  
  if(shortNames) {
    nm = colnames(geno)
    shortnm = ifelse(nchar(nm) > 7, paste0(substr(nm, 1, 4), "..."), nm)
    colnames(geno) = shortnm
  }
  
  # Allele labels: Sort and truncate
  alleleText = vapply(attrs, function(a) {
    labs = a$alleles
    nums = suppressWarnings(as.numeric(labs))
    labs = if(all(!is.na(nums))) labs[order(nums)] else sort.int(labs)
    shortAlleles(labs)
  }, "")

  # Mutation model names
  mods = vapply(attrs, \(a)
    if(is.null(a$mutmod)) "-" else getParams(a$mutmod, "model", format = 3)$model, "")

  # Main table
  x = data.frame(Marker = markers,
                 N = vapply(attrs, function(a) length(a$alleles), integer(1)),
                 Alleles = alleleText,
                 geno,
                 Mut = mods,
                 LR = formatLR(lr[markers]),
                 check.names = FALSE)

  if(!is.null(lrNoMut))
    x[["LR*"]] = formatLR(unname(lrNoMut[markers]))

  # Add signed percentage change to the masked table
  if(!is.null(referenceLR)) {
    deviation = rep(NA_character_, length(lr))
    ok = is.finite(lr) & is.finite(referenceLR) & referenceLR != 0
    deviation[ok] = sprintf("%+.1f%%", 100 * (lr[ok] / referenceLR[ok] - 1))
    x[[ncol(x) + 1]] = deviation
    names(x)[ncol(x)] = ""
  }

  DT::datatable(
    x, rownames = FALSE, class = "compact stripe hover nowrap",
    options = list(
      dom = "t", paging = FALSE, scrollX = TRUE, ordering = FALSE,
      scrollY = "330px", scrollCollapse = TRUE
    ),
    callback = DT::JS("markerDblclick(table);")
  )
}

plotAllPeds = function(peds) {
  npeds = length(peds)
  cex = if(npeds == 1) 1.1 else 1.3
  plotPedList(peds, hatched = typedMembers, cex = cex)
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
      scrollY = "60vh", scrollCollapse = TRUE,
      orderCellsTop = FALSE,
      columnDefs = list(list(
        targets = 0:3, className = "dt-center",
        render = DT::JS("sortBlanksLast")
      ))
    )
  )
}
