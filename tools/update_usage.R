# Regenerate public function usage sections from the package root.
# Rd descriptions are maintained separately; this updates usage only.
functions <- list()
for (path in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) {
  for (expr in parse(path)) {
    if (!is.call(expr) || length(expr) != 3L) next
    if (!as.character(expr[[1L]]) %in% c("<-", "=")) next
    if (!is.symbol(expr[[2L]])) next
    rhs <- expr[[3L]]
    if (!is.call(rhs) || !identical(rhs[[1L]], as.name("function"))) next
    functions[[as.character(expr[[2L]])]] <- eval(rhs, envir = baseenv())
  }
}

for (name in names(functions)) {
  if (length(formals(functions[[name]])) == 0L) next
  path <- file.path("man", paste0(name, ".Rd"))
  if (!file.exists(path)) next
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr("\\usage{", text, fixed = TRUE)[1L]
  end <- regexpr("\\arguments{", text, fixed = TRUE)[1L]
  if (start < 1L || end <= start) stop("Cannot locate usage in ", path)

  signature <- deparse(args(functions[[name]]), width.cutoff = 80L)
  if (!identical(tail(signature, 1L), "NULL")) stop("Unexpected signature: ", name)
  signature <- head(signature, -1L)
  signature[1L] <- sub("^function", name, signature[1L])
  signature <- paste(signature, collapse = "\n")
  for (character in c("\\", "{", "}", "%")) {
    signature <- gsub(character, paste0("\\", character), signature, fixed = TRUE)
  }
  updated <- paste0(substr(text, 1L, start - 1L), "\\usage{\n", signature,
                    "\n}\n", substring(text, end))
  # Validate before replacing a help file.
  connection <- textConnection(updated)
  tryCatch(
    withCallingHandlers(tools::parse_Rd(connection),
      warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    finally = close(connection)
  )
  writeLines(updated, path)
}
