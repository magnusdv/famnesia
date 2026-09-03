     
renderMarkerTable = function(peds, locusAttributes, lr, digits = 2, 
                             lrNoMut = NULL,
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
                 LR = formatLR(lr[markers], digits),
                 check.names = FALSE)

  if(!is.null(lrNoMut))
    x[["LR*"]] = formatLR(unname(lrNoMut[markers]), digits)

  # Add signed percentage change to the masked table
  if(!is.null(referenceLR)) {
    deviation = rep(NA_character_, length(lr))
    ok = is.finite(lr) & is.finite(referenceLR) & referenceLR != 0
    deviation[ok] = sprintf("%+.1f%%", 100 * (lr[ok] / referenceLR[ok] - 1))
    x[[ncol(x) + 1]] = deviation
    names(x)[ncol(x)] = ""
  }

  DT::datatable(
    x, rownames = FALSE, selection = "none", 
    class = "compact stripe hover nowrap",
    options = list(
      dom = "t", paging = FALSE, scrollX = TRUE, ordering = FALSE,
      scrollY = "300px", scrollCollapse = TRUE
    ),
    callback = DT::JS("markerDblclick(table);")
  ) |> 
    DT::formatStyle(1:ncol(x), lineHeight = "96%") 
}


# Marker-wise LRs
markerLR = function(peds, theta = 0) {
  peds12 = head(peds, 2)
  lik = do.call(cbind, lapply(peds12, likelihood, theta = theta))
  lr = if(ncol(lik) < 2) rep(NA_real_, nrow(lik)) else lik[, 1] / lik[, 2]
  setnames(lr, name(peds[[1]]))
}

