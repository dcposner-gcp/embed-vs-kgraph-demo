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
# PPM pinned to the DATED snapshot renv.lock was written against (2026-09-03):
# every locked version is current there, so all three OSes get prebuilt
# binaries. PPM's "latest" channel keeps binaries only for the newest release
# of each package, so as CRAN moves on a "latest" URL forces source builds of
# superseded pins (fatal without a compiler toolchain). Bump the date and
# re-snapshot together (same date in .github/workflows/smoke.yml). CRAN stays
# as a FALLBACK for anything the snapshot lacks.
# Linux needs the distro-specific __linux__ path (noble = the rocker Ubuntu
# 24.04 base); macOS/Windows use PPM's platform endpoint, which serves native
# binaries (the __linux__ path is Linux-only and would yield no Mac/Win builds).
if (Sys.info()[["sysname"]] == "Linux") {
  options(repos = c(
    PPM  = "https://packagemanager.posit.co/cran/__linux__/noble/2026-09-03",
    CRAN = "https://cloud.r-project.org"
  ))
} else {
  options(repos = c(
    PPM  = "https://packagemanager.posit.co/cran/2026-09-03",
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
