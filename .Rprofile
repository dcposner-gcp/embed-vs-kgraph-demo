# ============================================================================
# Project .Rprofile
# Sourced at every R session start (including non-interactive builds and
# RStudio launches). Keep this lean — heavy work belongs in scripts.
# ============================================================================

# --- renv activation -------------------------------------------------------
# Gated so the file is safe before renv has been initialized (fresh clone).
local({
  activate <- file.path("renv", "activate.R")
  if (file.exists(activate)) source(activate)
})

# --- Repos ------------------------------------------------------------------
# PPM (binaries) first, CRAN as a FALLBACK so packages PPM hasn't mirrored yet
# (it lags CRAN by a day or two) still resolve to the exact lockfile version.
# Linux needs the distro-specific __linux__ path (noble = the rocker Ubuntu
# 24.04 base); macOS/Windows use PPM's platform endpoint, which serves native
# binaries (the __linux__ path is Linux-only and would yield no Mac/Win builds).
if (Sys.info()[["sysname"]] == "Linux") {
  options(repos = c(
    PPM  = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
    CRAN = "https://cloud.r-project.org"
  ))
} else {
  options(repos = c(
    PPM  = "https://packagemanager.posit.co/cran/latest",
    CRAN = "https://cloud.r-project.org"
  ))
}

# --- pak engine for renv ---------------------------------------------------
# Parallelized installs via pak. pkgType="binary" is correct on macOS/Windows
# (CRAN ships binaries) but BREAKS on Linux: PPM has no arm64 binaries, so it
# serves source packages, and forcing "binary" makes pak reject them
# ("not a valid binary"). Only set it off-Linux; Linux falls back to source.
options(renv.config.pak.enabled = TRUE)
if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
  options(pkgType = "binary", install.packages.check.source = "no")
}

# --- Console / interactive niceties ----------------------------------------
options(
  digits          = 6,
  scipen          = 10,            # prefer 0.0001 over 1e-04
  width           = 120,
  warn            = 1,             # warnings as they happen, not at end
  stringsAsFactors = FALSE,        # belt-and-suspenders even on R >= 4.0
  datatable.print.class = TRUE,
  datatable.print.keys  = TRUE,
  datatable.print.trunc.cols = TRUE
)
