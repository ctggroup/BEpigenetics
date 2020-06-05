mytopGo <- function (results, ontology = c("BP", "CC", "MF"), sort = NULL, 
          number = 20L, truncate.term = NULL) 
{
  if (!is.data.frame(results)) 
    stop("results should be a data.frame.")
  ontology <- match.arg(unique(ontology), c("BP", "CC", "MF"), 
                        several.ok = TRUE)
  if (length(ontology) < 3L) {
    sel <- results$Ont %in% ontology
    results <- results[sel, ]
  }
  dimres <- dim(results)
  if (!is.numeric(number)) 
    stop("number should be a positive integer")
  if (number > dimres[1L]) 
    number <- dimres[1]
  if (number < 1L) 
    return(results[integer(0), ])
  nsets <- (dimres[2L] - 3L)%/%2L
  if (nsets < 1L) 
    stop("results has wrong number of columns")
  setnames <- colnames(results)[4L:(3L + nsets)]
  if (is.null(sort)) {
    isort <- 1L:nsets
  }
  else {
    sort <- as.character(sort)
    isort <- which(tolower(setnames) %in% tolower(sort))
    if (!length(isort)) 
      stop("sort column not found in results")
  }
  P.col <- 3L + nsets + isort
  if (length(P.col) == 1L) {
    P <- results[, P.col]
  }
  else {
    P <- do.call("pmin", as.data.frame(results[, P.col, 
                                               drop = FALSE]))
  }
  tmp <- results$Term
  o <- order(P, results$N)
  tab <- results[o[1L:number], , drop = FALSE]
  if (!is.null(truncate.term)) {
    truncate.term <- as.integer(truncate.term[1])
    truncate.term <- max(truncate.term, 5L)
    truncate.term <- min(truncate.term, 1000L)
    tm2 <- truncate.term - 3L
    i <- (nchar(tab$Term) > tm2)
    tab$Term[i] <- paste0(substring(tab$Term[i], 1L, tm2), 
                          "...")
  }
  tab
}