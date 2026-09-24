#' Look up reference numbers in the manuscript's own bibliography.
#'
#' Figure S6 prints citation numbers inside the plot. Hard-coding them means a
#' renumbered bibliography silently invalidates the figure, which has happened
#' twice. This reads the numbers out of the manuscript source instead, matching
#' each reference by first-author surname and year, and stops if a reference is
#' missing or matches more than one entry.
#'
#' @param keys Named list; each element is c(surname, year).
#' @param path Manuscript source containing a "## References" section.
#' @return Named integer vector of reference numbers, in the order of `keys`.
#' @examples
#'   n <- reference_numbers(list(korva = c("Korva", "2013")))
#'   sprintf("(Korva 2013 [%d])", n[["korva"]])
reference_numbers <- function(keys,
                              path = file.path("manuscript", "Royal", "_src",
                                               "manuscript.md.in")) {
  if (!file.exists(path)) {
    stop("cannot find the manuscript source at ", path,
         "; run this from the project root")
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  start <- grep("^## References", lines)
  if (length(start) != 1) {
    stop("expected exactly one '## References' heading, found ", length(start))
  }
  later <- grep("^## ", lines)
  later <- later[later > start]
  stop_at <- if (length(later) > 0) later[1] - 1 else length(lines)

  entries <- grep("^[0-9]+[.] ", lines[(start + 1):stop_at], value = TRUE)
  if (length(entries) == 0) stop("no numbered references found in ", path)

  out <- integer(0)
  for (nm in names(keys)) {
    surname <- keys[[nm]][1]
    year <- keys[[nm]][2]
    hit <- entries[grepl(surname, entries, fixed = TRUE) &
                     grepl(year, entries, fixed = TRUE)]
    if (length(hit) != 1) {
      stop(sprintf("expected one reference matching '%s %s', found %d",
                   surname, year, length(hit)))
    }
    out[nm] <- as.integer(sub("^([0-9]+)[.].*$", "\\1", hit))
  }
  out
}
