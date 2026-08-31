#!/usr/bin/env Rscript
# scripts/simulate_data.R -- generate a complete simulated `data/` tree.
#
# Writes every input `_targets.R` reads, from one latent geometry with a known
# signal (functions in R/simulate.R), so the unchanged pipeline runs end to end
# with no private data:
#   <out>/dim/Phecode_map_v1_2_icd9_icd10cm.csv     public Phecode 1.2 (PheWAS mirror)
#   <out>/dim/phecode_definitions1.2.csv            public Phecode 1.2
#   <out>/dim/icd_athena_2025/CONCEPT.csv           SYNTHETIC ICD typing dictionary
#   <out>/validation/dict-disease-disease-wiki-bin-eval.Rdata   simulated gold graph
#   <out>/prediction/<query>_Diagnosis.csv          simulated retriever-shaped retrievals
#   <out>/spec/{nodes,gold_edges,queries,predictions,embeddings}.csv  SPEC section-1 tables (docs/SPEC-inputs.md)
#   <out>/sim_manifest.json, <out>/README.md         parameters, seed, truth, provenance
#
# Usage (from the project root):
#   Rscript scripts/simulate_data.R                       # preset phecode -> data_sim/phecode
#   Rscript scripts/simulate_data.R --preset phecode-null # snr = 0: AUROC should be ~0.5
#   Rscript scripts/simulate_data.R --preset phecode-oracle --out data_sim/oracle
#   Rscript scripts/simulate_data.R --set snr=4 --set lambda=0.9 --seed 7
#   Rscript scripts/simulate_data.R --out data --force    # fresh clone: write straight to data/
# Then run the pipeline: scripts/run_pipeline_sim.sh <out>   (or tar_make() if <out> is data/)
#
# Refuses to overwrite an existing prediction/ directory without --force, so a
# real data/ tree is never clobbered by accident.

suppressPackageStartupMessages({
  library(data.table)
  library(cli)
})

# ---- locate project root + source the simulator ----------------------------
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- if (length(script_path)) normalizePath(file.path(dirname(script_path), "..")) else getwd()
setwd(root)
source(file.path(root, "R", "simulate.R"))

# ---- presets -----------------------------------------------------------------
Presets <- function() {
  base <- list(
    # latent comorbidity view V (generates the gold graph)
    d = 32, sd_category = 1, sd_root = 1, sd_mid = 0.7, sd_leaf = 0.5,
    # hierarchy / lexical view H (what a description embedder sees)
    h_sd_category = 0.8, h_sd_root = 1, h_sd_mid = 0.7, h_sd_leaf = 0.15,
    # gold graph (closure = 0.1 gives ~40% hierarchy compliance under the pipeline's
    # CheckPhecodeHierarchy() definition; the real KESER graph shows ~34%)
    edge_rule = "logistic", mean_degree = 12, beta = 6, gamma = 2.5, sd_pop = 0.8,
    closure = 0.1, parent_dropout = 0.3, noise_frac = 0.05, min_degree = 5,
    zero_positive = "496.21", force_pairs = "411.2:411.4,250.2:411.4",
    # edge_rule = "loglinear" only (Xu et al. generative model; see R/simulate.R).
    # Density is emergent: 5k patients x 100 events at kappa 1.5 -> mean degree ~14
    # (the real KESER baseline's query degrees are ~15); 20k -> ~50.
    n_patients = 5000, seq_len = 100, kappa = 1.5, fdr = 0.05, min_count = 5,
    # retriever: U = lambda*V + (1-lambda)*H + noise/snr; top-k phecodes per query
    lambda = 0.7, snr = 1.5, k_mean = 1200, k_sdlog = 0.3,
    # ICD layer quirks
    lambda_extra = 0.8, sd_icd = 0.01, p_list = 0.04, p_range = 0.01,
    p_ambiguous = 0.02, p_unmapped = 0.05, p_desc_perturb = 0.3, ambiguity_frac = 0.5,
    icd_unique_only = FALSE
  )
  list(
    phecode = base,
    `phecode-null` = modifyList(base, list(snr = 0)),
    # gold graph from simulated patient sequences (Xu et al. 2026); ~1 min at 5k patients
    `phecode-xu` = modifyList(base, list(edge_rule = "loglinear")),
    `phecode-oracle` = modifyList(base, list(
      edge_rule = "threshold", gamma = 0, sd_pop = 0, closure = 0, parent_dropout = 0,
      noise_frac = 0, lambda = 1, snr = Inf, k_mean = Inf, sd_icd = 0, p_range = 0,
      ambiguity_frac = 0, icd_unique_only = TRUE, force_pairs = ""))
  )
}

# The query set is fixed by LoadEmbedDx() in R/functions.R (file stem ->
# Phecode 1.2). COVID / Long_COVID have no 1.2 code: their retrieval is anchored
# on a stand-in node and their rank-1 self code is an unmapped ICD-10 code.
QueryTable <- function() data.table(
  query_id  = c("CAD", "COPD", "COVID", "Depression", "Hearing_Loss", "HF",
                "Long_COVID", "MI", "MS", "RA", "T2DM"),
  label     = c("Coronary atherosclerosis", "COPD", "COVID-19", "Major depressive disorder",
                "Hearing loss", "Heart failure", "Long COVID", "Myocardial infarction",
                "Multiple sclerosis", "Rheumatoid arthritis", "Type 2 diabetes"),
  node_id   = c("411.4", "496.21", NA, "296.2", "389", "428", NA, "411.2", "335", "714.1", "250.2"),
  anchor    = c(NA, NA, "480", NA, NA, NA, "798", NA, NA, NA, NA),
  self_code = c(NA, NA, "U07.1", NA, NA, NA, "U09.9", NA, NA, NA, NA)
)

# ---- CLI -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(preset = "phecode", out = NULL, seed = 20260829L,
            cache = "data_sim/_source", force = FALSE, set = list())
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  nxt <- function() { i <<- i + 1L; if (i > length(args)) stop("missing value after ", a); args[i] }
  switch(a,
    "--preset" = opt$preset <- nxt(),
    "--out"    = opt$out    <- nxt(),
    "--seed"   = opt$seed   <- as.integer(nxt()),
    "--cache"  = opt$cache  <- nxt(),
    "--force"  = opt$force  <- TRUE,
    "--set"    = {
      kv <- strsplit(nxt(), "=", fixed = TRUE)[[1L]]
      if (length(kv) != 2L) stop("--set expects key=value, got: ", args[i])
      num <- suppressWarnings(as.numeric(kv[2L]))
      opt$set[[kv[1L]]] <- if (is.na(num)) kv[2L] else num
    },
    "--help" = , "-h" = { cat(readLines(script_path, n = 25L)[-1L], sep = "\n"); quit(status = 0) },
    stop("unknown argument: ", a)
  )
  i <- i + 1L
}
presets <- Presets()
if (!opt$preset %in% names(presets))
  stop("unknown preset '", opt$preset, "'; choose one of: ", paste(names(presets), collapse = ", "))
p <- modifyList(presets[[opt$preset]], opt$set)
unknown <- setdiff(names(opt$set), names(presets$phecode))
if (length(unknown)) stop("unknown --set key(s): ", paste(unknown, collapse = ", "))
out <- if (is.null(opt$out)) file.path("data_sim", opt$preset) else opt$out

# ---- safety: never clobber an existing data tree silently ---------------------------
pred_dir <- file.path(out, "prediction")
if (dir.exists(pred_dir) && length(list.files(pred_dir)) && !opt$force)
  stop(sprintf("%s already has prediction files; pass --force to overwrite (NEVER on a real data/ tree).", out))
if (normalizePath(out, mustWork = FALSE) == normalizePath("data", mustWork = FALSE) && !opt$force)
  stop("refusing to write into data/ without --force (intended for a fresh clone only).")

cli_h1("simulate_data.R -- preset {opt$preset}, seed {opt$seed} -> {out}")

# ---- 1. vocabulary (public Phecode 1.2) ---------------------------------------------
src <- FetchPhecodeSource(opt$cache)
nodes <- BuildSimVocabulary(src)
queries <- QueryTable()
missing_q <- queries[!is.na(node_id) & !(node_id %chin% nodes$node_id), node_id]
if (length(missing_q)) stop("query phecode(s) not in vocabulary: ", paste(missing_q, collapse = ", "))
queries[is.na(node_id) & !(anchor %chin% nodes$node_id), anchor := nodes$node_id[1L]]
# self code for coded queries: the shortest ICD-9-CM code of the query phecode
self9 <- src$map[vocabulary_id == "ICD9CM"][order(nchar(code), code)][, .SD[1L], by = phecode]
queries[self9, on = .(node_id = phecode), self_code := fifelse(is.na(self_code), i.code, self_code)]
cli_alert_success("{nrow(nodes)} phecodes, {nrow(src$map)} ICD codes, {uniqueN(nodes$category)} categories")

# ---- 2. latent geometry --------------------------------------------------------------
V <- SimulateLatent(nodes, d = p$d, sd_category = p$sd_category, sd_root = p$sd_root,
                    sd_mid = p$sd_mid, sd_leaf = p$sd_leaf, seed = opt$seed + 1L)
H <- SimulateLatent(nodes, d = p$d, sd_category = p$h_sd_category, sd_root = p$h_sd_root,
                    sd_mid = p$h_sd_mid, sd_leaf = p$h_sd_leaf, seed = opt$seed + 2L)

# ---- 3. gold graph -------------------------------------------------------------------
query_nodes <- queries[!is.na(node_id), node_id]
zero_pos <- if (nzchar(p$zero_positive)) strsplit(p$zero_positive, ",", fixed = TRUE)[[1L]] else character()
force_pairs <- if (nzchar(p$force_pairs)) {
  lapply(strsplit(p$force_pairs, ",", fixed = TRUE)[[1L]], function(s) strsplit(s, ":", fixed = TRUE)[[1L]])
} else {
  list()
}
g <- SimulateGoldEdges(nodes, V, edge_rule = p$edge_rule, mean_degree = p$mean_degree,
                       beta = p$beta, gamma = p$gamma, sd_pop = p$sd_pop, closure = p$closure,
                       parent_dropout = p$parent_dropout, noise_frac = p$noise_frac,
                       protect = query_nodes, seed = opt$seed + 3L,
                       n_patients = as.integer(p$n_patients), seq_len = as.integer(p$seq_len),
                       kappa = p$kappa, fdr = p$fdr, min_count = as.integer(p$min_count))
edges <- ApplyQueryConstraints(g$edges, V, query_nodes, zero_positive = zero_pos,
                               force_pairs = force_pairs, min_degree = as.integer(p$min_degree))
deg <- rbind(edges[, .(node = u)], edges[, .(node = v)])[, .N, by = node]
cli_alert_success(paste0("gold graph: {nrow(edges)} edges (",
  paste(sprintf("%s=%d", names(table(edges$src)), as.integer(table(edges$src))), collapse = ", "),
  "), {uniqueN(c(edges$u, edges$v))} nodes, {length(g$dropped_parents)} parents dropped"))

# ---- 4. retrieval ----------------------------------------------------------------------
U <- SimulateRetrieverEmbedding(V, H, lambda = p$lambda, snr = p$snr, seed = opt$seed + 4L)
qvec <- queries[, fifelse(is.na(node_id), anchor, node_id)]
Q <- U[qvec, , drop = FALSE]
rownames(Q) <- queries$query_id
self_node <- setNames(as.list(queries$node_id), queries$query_id)
pred <- ScorePhecodeRetrieval(U, Q, self_node, k_mean = p$k_mean, k_sdlog = p$k_sdlog,
                              seed = opt$seed + 5L)

# ---- 5. ICD layer ----------------------------------------------------------------------
dict <- BuildSimIcdDictionary(src$map, ambiguity_frac = p$ambiguity_frac, seed = opt$seed + 6L)
icd_rows <- SimulateIcdLayer(pred, queries, src, dict$ambiguous_codes,
                             lambda_extra = p$lambda_extra, sd_icd = p$sd_icd, p_list = p$p_list,
                             p_range = p$p_range, p_ambiguous = p$p_ambiguous,
                             p_unmapped = p$p_unmapped, p_desc_perturb = p$p_desc_perturb,
                             icd_unique_only = isTRUE(as.logical(p$icd_unique_only)),
                             seed = opt$seed + 7L)

# ---- 6. write --------------------------------------------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
files <- c(WriteSimPhecodeFiles(src, out), WriteSimConceptCsv(dict$dict, out),
           WriteSimGoldRdata(edges, out), WriteSimPredictions(icd_rows, out))
spec_dir <- WriteSimSpecTables(nodes, edges, queries, pred, out, U = U)

# ---- 7. truth + manifest (summaries only; no rows) -----------------------------------------
per_query <- queries[, {
  q <- query_id; nid <- node_id
  pos <- if (is.na(nid)) character() else edges[u == nid | v == nid, fifelse(u == nid, v, u)]
  got <- pred[query_id == q, candidate_id]
  .(node_id = nid, n_pos = length(pos), k_phecodes = length(got),
    recall_at_k = if (length(pos)) round(mean(pos %chin% got), 4) else NA_real_,
    n_icd_rows = icd_rows[query_id == q, .N])
}, by = query_id]
compliance <- tryCatch({
  # pipeline definition of hierarchy closure compliance, if the functions are available
  source(file.path(root, "R", "functions.R"), local = TRUE)
  e <- new.env(); load(file.path(out, "validation", "dict-disease-disease-wiki-bin-eval.Rdata"), envir = e)
  qc <- CheckPhecodeHierarchy(as.data.table(e$dict.wiki.keep))
  qc[metric == "compliance_pct", value]
}, error = function(err) NA_real_, warning = function(w) suppressWarnings({
  source(file.path(root, "R", "functions.R"), local = TRUE)
  e <- new.env(); load(file.path(out, "validation", "dict-disease-disease-wiki-bin-eval.Rdata"), envir = e)
  CheckPhecodeHierarchy(as.data.table(e$dict.wiki.keep))[metric == "compliance_pct", value]
}))
manifest <- list(
  generator = "scripts/simulate_data.R (R/simulate.R)",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  git_sha = tryCatch(trimws(system("git rev-parse HEAD", intern = TRUE)), error = function(e) NA),
  preset = opt$preset, seed = opt$seed, params = p,
  source = src$spec[c("map", "info", "citation")],
  vocabulary = list(n_nodes = nrow(nodes), n_icd_codes = nrow(src$map),
                    n_categories = uniqueN(nodes$category)),
  gold = list(n_edges = nrow(edges), edges_by_source = as.list(table(edges$src)),
              n_nodes_with_edges = uniqueN(c(edges$u, edges$v)),
              mean_degree = round(mean(deg$N), 3),
              n_parents_dropped = length(g$dropped_parents),
              hierarchy_compliance_pct = compliance,
              loglinear = if (!is.null(g$loglinear)) g$loglinear$stats else NULL),
  retrieval = list(icd_rows_by_kind = as.list(table(icd_rows$kind))),
  per_query = per_query,
  files = c(files, spec_dir)
)
jsonlite::write_json(manifest, file.path(out, "sim_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

writeLines(c(
  "# SIMULATED DATA -- nothing in this directory is real",
  "",
  sprintf("Generated by `scripts/simulate_data.R --preset %s --seed %d` (git %s).",
          opt$preset, opt$seed, substr(manifest$git_sha, 1, 7)),
  "",
  "| file | what it stands in for | provenance |",
  "|---|---|---|",
  "| `dim/Phecode_map_v1_2_icd9_icd10cm.csv`, `dim/phecode_definitions1.2.csv` | Phecode 1.2 map + definitions | **public**: PheWAS R package mirror (Denny lab); `ICDString` = phecode description |",
  "| `dim/icd_athena_2025/CONCEPT.csv` | OHDSI/Athena ICD vocabulary | **synthetic** from the public map (see README.txt there) |",
  "| `validation/dict-disease-disease-wiki-bin-eval.Rdata` | KESER dx-dx gold graph | **simulated** from latent vectors V |",
  "| `prediction/<query>_Diagnosis.csv` | embedding-model ranked ICD retrievals | **simulated** from U = lambda*V + (1-lambda)*H + noise/snr |",
  "| `spec/*.csv` | SPEC section-1 tables (docs/SPEC-inputs.md): nodes, gold_edges, queries, predictions, embeddings = the retriever vectors U | same simulation, phecode level |",
  "| `sim_manifest.json` | parameters, seed, per-query truth (n_pos, recall@k) | -- |",
  "",
  "Run the pipeline on it: `scripts/run_pipeline_sim.sh <this dir>` (isolated run",
  "dir) or, on a fresh clone, generate with `--out data --force` and `targets::tar_make()`."
), file.path(out, "README.md"))

cli_h2("summary")
print(per_query, row.names = FALSE)
cat(sprintf("\nhierarchy compliance (pipeline definition): %s%%\n", compliance))
cat("ICD rows by kind:", paste(sprintf("%s=%d", names(table(icd_rows$kind)), as.integer(table(icd_rows$kind))), collapse = ", "), "\n")
cli_alert_success("wrote {length(files)} pipeline input files + spec/ + sim_manifest.json to {out}")
