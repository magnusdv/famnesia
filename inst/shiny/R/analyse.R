analyseFam = function(x) {
  base = c("famnames", "ids", "markernames", "shuffle", "sex")
  mods = unlist(lapply(x$mutpars, \(p) p$model %||% character()), use.names = FALSE)
  stepwise = "stepwise" %in% mods
  hasMut = length(mods) > 0

  info = list(
    filename = x$filename,
    famnames = names(x$peds),
    individuals = length(typedMembers(x$peds[[1]])),
    database = sprintf("%s (%d markers)", x$params$dbName %||% "unnamed", length(x$attrs)),
    mutmodel = if(!hasMut) "None" else paste(unique(mods), collapse = "/")
  )

  if(!hasMut) {
    suggested = list(
      options = c(base, "lump"),
      alleles = "strong",
      freqs = "tweak",
      mutmodels = "original"
    )
    preserve = suggested
    preserve$freqs = "original"
    effect = "near-exact"
  }
  else if(!stepwise) {
    suggested = list(
      options = c(base),
      alleles = "strong",
      freqs = "tweak",
      mutmodels = "original"
    )
    preserve = list(
      options = base,
      alleles = "strong",
      freqs = "original",
      mutmodels = "original"
    )
    effect = "near-exact"
  }
  else {
    suggested = list(
      options = c(base),
      alleles = "strong",
      freqs = "tweak",
      mutmodels = "simplify"
    )
    preserve = list(
      options = base,
      alleles = "constrained",
      freqs = "original",
      mutmodels = "original"
    )
    effect = "non-exact"
  }

  # Maximal masking (always the same)
  max = list(
    options = c("famnames", "ids", "markernames", "shuffle", "lump", "sex"),
    alleles = "strong", freqs = "tweak", mutmodels = "disable"
  )
  
  list(info = info, max = max, suggested = suggested,
       preserve = preserve, effect = effect)
}

effectSymbol = function(effect) {
  switch(effect, "non-exact" = "≠", "near-exact" = "≈", "exact" = "=", "?")
}

maskingText = function(x) {
  txt = c(
    if("lump" %in% x$options) "lump alleles",
    switch(x$alleles,
           strong = "randomise alleles",
           constrained = "constrained alleles",
           NULL),
    switch(x$freqs,
           original = "original freqs",
           round = "round freqs",
           tweak = "tweak freqs"),
    switch(x$mutmodels,
           original = "original mut models",
           simplify = "simplify mut models",
           disable = "disable mut models")
  )
  paste(txt, collapse = ", ")
}

