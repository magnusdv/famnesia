maskData = function(x, famnames = FALSE, ids = FALSE,
                    markernames = FALSE, lump = FALSE, sex = FALSE,
                    alleles = c("original", "constrained","strong"),
                    freqs = c("original", "tweak", "round"),
                    mutmodels = c("original", "simplify", "disable")) {
  
  alleles = match.arg(alleles)
  freqs = match.arg(freqs)
  mutmodels = match.arg(mutmodels)

  # Starting point: Original peds without mut models
  peds = lapply(x$peds, setMutmod, model = NULL)

  # Mask: Randomise sex of suitable individuals
  if(sex)
    peds = randomizeSex(peds)

  # Mask: Allele lumping
  if(lump)
    peds = lapply(peds, \(p) lumpAlleles(p, observed = x$observed, force = TRUE))

  # Extract frequency database
  db = getFreqDatabase(peds[[1]], format = "list")
  
  # Mask: Tweak or round frequencies
  db = lapply(db, maskFreqs, method = freqs)
  
  # Mask: Allele labels
  oldAlleles = lapply(db, names)
  db = lapply(db, maskAlleleLabs, method = alleles)
  alleleMap = Map(\(old, new) setnames(new, old),
                  oldAlleles, lapply(db, names))
  mnames = names(db) |> setnames()
  
  # Mask: Marker names
  newMnames = if(markernames) sample(paste0("M", seq_along(mnames))) else mnames
  names(newMnames) = mnames
  
  # Rebuild mut models
  mutmods = lapply(mnames, \(m)
    restoreMutmod(db[[m]], x$mutpars[[m]], method = mutmodels))

  # Rebuild marker attributes
  attrs = lapply(mnames, \(m)
    list(name = newMnames[[m]], alleles = names(db[[m]]),
         afreq = unname(db[[m]]), mutmod = mutmods[[m]]))
  names(attrs) = unname(newMnames)

  # Apply the new labels to the integer genotype data
  peds = lapply(peds, setMaskedAttribs, attrs = attrs)

  # Mask: ID labels
  idList = lapply(peds, labels)
  oldids = unique.default(unlist(idList, use.names = FALSE))
  newids = if(ids) as.character(seq_along(oldids)) else oldids
  names(newids) = oldids
  if(ids)
    peds = lapply(seq_along(peds), \(i) relabel(peds[[i]], new = newids[idList[[i]]]))

  # Mask: Family names
  if(famnames)
    names(peds) = paste0("F", seq_along(peds))
  
  lrNoMut = if(!is.null(x$lrNoMut))
    markerLR(lapply(peds, setMutmod, model = NULL), theta = 0)

  list(
    peds = peds,
    attrs = attrs,
    alleleMap = alleleMap,
    lr = markerLR(peds, theta = 0),
    lrNoMut = lrNoMut,
    params = maskFamParams(x$params, newids, newMnames)
  )
}


setMaskedAttribs = function(x, attrs) {
  if(is.pedList(x))
    return(lapply(x, setMaskedAttribs, attrs = attrs))
  
  m = as.matrix(x, include.attrs = TRUE)
  a = attributes(m)
  a$markerattr = Map(modifyList, a$markerattr, attrs)
  restorePed(m, attrs = a)
}


maskFreqs = function(v, method) {
  if(method == "original")
    return(v)
  
  switch(method,
    tweak = {
      vv = v * runif(length(v), 0.98, 1.02)
      vv / sum(vv)
   },
    round = {
      v[v < 0.0001] = 0.0001
      vv = round(v, 4)
      j = which.max(vv)
      vv[j] = 1 - sum(vv[-j])
      round(vv, 4)
    }
  )
}

maskAlleleLabs = function(v, method) {
  if(method == "original")
    return(v)
  
  als = names(v)
  newals = switch(method, 
    strong = as.character(sample.int(length(als))),
    constrained = constrainedLabels(als)
  )
  setnames(v, newals)
}

restoreMutmod = function(afreq, pars, method) {
  if(method == "disable" || is.null(pars))
    return(NULL)
  
  if(method == "simplify")
    pars$model = "equal"
  
  do.call(pedmut::mutationModel, c(pars, list(afreq = afreq)))
}

randomizeSex = function(peds) {
  # Keep IDs that are leaves in every pedigree
  ids = Reduce(intersect, lapply(peds, leaves))
  newsex = sample(1:2, length(ids), replace = TRUE) |>
    setnames(ids)
  lapply(peds, setSex, sex = newsex)
}

maskFamParams = function(params, idLabs = NULL, markerNames = NULL) {
  p = params
  p$dbName = "<masked>"

  # ID labels in dropout parameters
  if(!is.null(p$dropoutConsider) && !is.null(idLabs))
    names(p$dropoutConsider) = unname(idLabs[names(p$dropoutConsider)])

  # Marker names in various parameters
  if(!is.null(markerNames)) {
    oldnames = names(markerNames)
    newnames = unname(markerNames)
    markerParams = intersect(c("dbSize", "dropoutValue", "maf"), names(p))
    p[markerParams] = lapply(p[markerParams], function(v)
      if(is.null(v)) NULL else setnames(v[oldnames], newnames))
  }
  
  p
}


constrainedLabels = function(alleles) {
  old = as.character(alleles)
  x = suppressWarnings(as.numeric(old))
  isNumeric = !is.na(x)

  if(!any(isNumeric))
    return(as.character(seq_along(old)))

  whole = floor(x[isNumeric])
  decimal = round(x[isNumeric] - whole, 10)
  groups = sort(unique(decimal[decimal != 0]))
  if(length(groups) > 9)
    stop("Constrained masking supports at most nine decimal groups", call. = FALSE)

  decimal[decimal != 0] = match(decimal[decimal != 0], groups) / 10
  new = numeric(length(old))
  new[isNumeric] = whole - min(whole) + 1 + decimal

  # Put nonnumeric labels such as the lump after the numeric ladder
  new[!isNumeric] = floor(max(new[isNumeric])) + seq_len(sum(!isNumeric))

  labels = format(round(new, 1), trim = TRUE, scientific = FALSE, nsmall = 1)
  labels = sub("\\.0$", "", labels)
  labels
}
