cat("Loading custom .Rprofile\n")

## Set CRAN mirror automatically
.First <- function() {
    options(
        repos = c(CRAN = "https://muug.ca/mirror/cran/")
	)
}

## Increase size of Rhistory file
Sys.setenv(R_HISTSIZE='100000')

## Set Library
# .libPaths(c('~/bin/R/x86_64-pc-linux-gnu-library/4.0', .libPaths()))
.libPaths(c('~/bin/R/x86_64-pc-linux-gnu-library/4.4-paga', .libPaths()))
cat("\n")
print(.libPaths())
cat("\n")
cat("Library paths modified\n")
