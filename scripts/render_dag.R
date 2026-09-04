#!/usr/bin/env Rscript
# =============================================================================
# render_dag.R  —  Draw the `_targets.R` dependency graph as a static figure.
# -----------------------------------------------------------------------------
# WHAT: docs/figures/pipeline_dag.svg — every data target in the pipeline as a
#   box, edges meaning "feeds into", coloured by pipeline stage (the numbered
#   stages in README.md). Targets that track files (raw inputs, written outputs
#   and QC tables) are drawn as file shapes. The GraphViz source goes to
#   output/dag/pipeline_dag.dot (gitignored intermediate).
#
# WHY a script, not a target: the figure documents the pipeline; nothing in the
#   pipeline depends on it. `targets::tar_network()` reads _targets.R only — no
#   built store, no simulated data — so this runs on a fresh clone. It colours
#   by stage rather than by build status (what `tar_visnetwork()` does) because
#   a README figure should describe the pipeline, not whether one machine's
#   store is current; on a fresh clone every node would be "outdated" anyway.
#
# RUN:  Rscript scripts/render_dag.R     (from the repo root)  |  make dag
# Needs the project's R library (renv::restore) and the GraphViz `dot` binary
# (`brew install graphviz` / `apt install graphviz`). Without `dot` the .dot
# source is still written and the SVG step is skipped with a hint.
# =============================================================================

if (!file.exists("_targets.R")) {
  stop("run from the repo root (no _targets.R in ", getwd(), ")", call. = FALSE)
}
suppressPackageStartupMessages(library(targets))

# --- 1. The graph: data targets only ------------------------------------------
net <- tar_network(targets_only = TRUE, reporter = "silent")
v   <- net$vertices
e   <- net$edges
# Which targets track files (format = "file"): drawn as file shapes below.
fmt <- tar_manifest(fields = c("name", "format"))
is_file <- v$name %in% fmt$name[fmt$format == "file"]

# --- 2. Stage of each target --------------------------------------------------
# One fill per role, the same palette as docs/figures/concept.dot (blue = the
# knowledge-graph answer key, pink = the embedding under test, grey = synthetic
# controls, gold = results, neutral = inputs). A target missing from this table
# is drawn white and reported, so a new target in _targets.R gets noticed.
stage_of <- c(
  # raw inputs and the dictionaries loaded from them
  phe_file = "input", icd_dict_file = "input", keser_dx_file = "input",
  predict_dir = "input", phecode_def_file = "input",
  phemap = "input", dict_icd = "input", phe_nodes = "input",
  # the embedding's retrievals, typed and mapped to phecodes, plus their QC
  embed_clinical = "retrieval", audit_results = "retrieval",
  write_audit_tsv = "retrieval", write_audit_xlsx = "retrieval",
  qc_range_codetypes = "retrieval", phecodes_to_validate = "retrieval",
  # the knowledge-graph gold standard
  dict_wiki_keep = "gold", keser_hierarchical = "gold",
  write_keser_hierarchy = "gold", dict_keser_full = "gold", dt_true = "gold",
  # matched synthetic controls
  dt_synthetic_controls = "controls",
  # scoring, metrics and the written outputs
  dt_predicted = "results", dt_eval = "results", dt_scored = "results",
  dt_auc = "results", dt_roc = "results", plot_auc = "results",
  query_phecode_colors = "results", query_phecode_labels = "results",
  out_core = "results", out_auc_by_query = "results",
  out_roc_by_query = "results", out_headline_auc = "results"
)
palette <- list(
  input     = list(fill = "#fafafa", line = "#8a8a8a", legend = "inputs and dictionaries"),
  retrieval = list(fill = "#f0d8e6", line = "#b0588f", legend = "embedding retrievals (under test)"),
  gold      = list(fill = "#d5e6f5", line = "#2b7bb9", legend = "knowledge graph (gold standard)"),
  controls  = list(fill = "#eeeeee", line = "#8a8a8a", legend = "synthetic controls"),
  results   = list(fill = "#fdf3e3", line = "#c9a24d", legend = "scoring, metrics, outputs"),
  unknown   = list(fill = "#ffffff", line = "#d32f2f", legend = "not in the stage table")
)
stage <- unname(stage_of[v$name])
stage[is.na(stage)] <- "unknown"
if (any(stage == "unknown")) {
  message("! targets missing from the stage table in scripts/render_dag.R: ",
          paste(v$name[stage == "unknown"], collapse = ", "))
}
stale <- setdiff(names(stage_of), v$name)
if (length(stale)) {
  message("! stage table lists targets that no longer exist: ", paste(stale, collapse = ", "))
}
pick <- function(field) vapply(stage, function(s) palette[[s]][[field]], character(1))

# --- 3. GraphViz source -------------------------------------------------------
esc <- function(x) gsub('"', '\\\\"', x)
node_lines <- sprintf(
  '  "%s" [label="%s", shape=%s, style="%s", fillcolor="%s", color="%s"];',
  esc(v$name), esc(v$name),
  ifelse(is_file, "note", "box"),
  ifelse(is_file, "filled", "filled,rounded"),
  pick("fill"), pick("line")
)
edge_lines <- if (nrow(e)) sprintf('  "%s" -> "%s";', esc(e$from), esc(e$to)) else character(0)
used <- names(palette)[names(palette) %in% unique(stage)]
legend_nodes <- sprintf(
  '    L_%s [label="%s", shape=box, style="filled,rounded", fillcolor="%s", color="%s"];',
  used, vapply(used, function(s) palette[[s]]$legend, character(1)),
  vapply(used, function(s) palette[[s]]$fill, character(1)),
  vapply(used, function(s) palette[[s]]$line, character(1))
)
legend <- c(
  '  subgraph cluster_legend {',
  '    label="legend (file shape = a tracked file)"; fontsize=9.5; color="#bbbbbb"; style="rounded";',
  '    node [fontsize=9];',
  legend_nodes,
  sprintf('    %s [style=invis];', paste(paste0("L_", used), collapse = " -> ")),
  '  }'
)
dot_src <- c(
  "// GENERATED by scripts/render_dag.R from _targets.R — edit the script, not this file.",
  "digraph pipeline_dag {",
  '  rankdir=LR;',
  # solid background: GitHub's dark theme renders a transparent SVG as dark text on dark
  '  graph [bgcolor="white", fontname="Helvetica,Arial,sans-serif", nodesep=0.25, ranksep=0.6];',
  '  node  [fontname="Helvetica,Arial,sans-serif", fontsize=10, penwidth=1];',
  '  edge  [color="#9e9e9e", arrowsize=0.7];',
  node_lines, edge_lines, legend, "}"
)

# --- 4. Write the source, render the SVG --------------------------------------
dot_dir  <- file.path("output", "dag")
dir.create(dot_dir, showWarnings = FALSE, recursive = TRUE)
dot_path <- file.path(dot_dir, "pipeline_dag.dot")
writeLines(dot_src, dot_path)
message(sprintf("graphviz source -> %s  (%d targets, %d edges)", dot_path, nrow(v), nrow(e)))

fig_dir <- file.path("docs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
svg_path <- file.path(fig_dir, "pipeline_dag.svg")
dot_bin  <- Sys.which("dot")
if (!nzchar(dot_bin)) {
  message("! graphviz `dot` not found: wrote the .dot source but skipped ", svg_path)
  message("  install it (brew install graphviz / apt install graphviz) and re-run.")
} else {
  status <- system2(dot_bin, c("-Tsvg", shQuote(dot_path), "-o", shQuote(svg_path)))
  if (status == 0L && file.exists(svg_path)) {
    message("static figure -> ", svg_path)
  } else {
    stop("`dot` failed (exit ", status, "); SVG not written", call. = FALSE)
  }
}
