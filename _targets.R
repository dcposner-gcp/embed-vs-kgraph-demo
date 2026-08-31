library(targets)
library(tarchetypes)

tar_source()
tar_option_set(
  packages = c("data.table", "stringi", "stringdist", "openxlsx", "ggplot2"),
  format = "rds"
)

# Single source of truth for evaluation output paths. Namespaced by retrieval
# task: "dx" = Diagnoses -> Diagnoses (Phecode). Add "output/rx" (a second
# OUT_DIR / pipeline branch) when the Diagnoses -> Medications task lands.
# Change this one line to repoint every output sink below.
OUT_DIR <- "output/dx"

list(
  # ./data/dim/ = dim and dictionary files
  tar_target(phe_file, "data/dim/Phecode_map_v1_2_icd9_icd10cm.csv", format = "file"),
  tar_target(icd_dict_file, "data/dim/icd_athena_2025/CONCEPT.csv", format = "file"),
  # ./data/validation/ = validation data for evaluating predictions
  tar_target(keser_dx_file,
    "data/validation/dict-disease-disease-wiki-bin-eval.Rdata", format = "file"),
  # ./data/prediction/ = prediction data
  # tracks all files in directory. any change will re-run entire pipeline
  tar_target(predict_dir, "data/prediction/", format = "file"),
  # load Denny phecode map 1.2
  tar_target(phemap, LoadPhecodeMap(phe_file)),
  # load phecode definitions
  tar_target(phecode_def_file, "data/dim/phecode_definitions1.2.csv", format = "file"),
  # load Denny phecode definitions
  tar_target(phe_nodes,
    {
      unique(fread(phecode_def_file,
            select = c("phecode" = "character", "category_number" = "integer")))
    }
  ),
  # ICD-9/10CM dictionary from Athena OHDSI
  tar_target(dict_icd, PopulateDictionary(icd_dict_file)),
  # load KESER dx-dx pairs (gold std)
  tar_target(dict_wiki_keep, LoadKeserDxDx(keser_dx_file)),
  # check KESER phecode hierarchy closure
  tar_target(keser_hierarchical, CheckPhecodeHierarchy(dict_wiki_keep)),
  # [human_review]: hierarchy QC table  
  tar_target(write_keser_hierarchy,
    {
      dir.create("scratch", showWarnings = FALSE, recursive = TRUE)
      out_path <- "scratch/keser_hierarchy_qc.tsv"
      fwrite(keser_hierarchical, out_path, sep = "\t")
      return(out_path)
    },
    format = "file"
  ),
  # load embedding retrievals, split lists & ranges, get ICD type, map to phecodes
  tar_target(embed_clinical, ProcessClinicalData(predict_dir, dict_icd, phemap)),
  # audit ambiguous ICD codes
  tar_target(audit_results, GenerateAmbiguityAudit(embed_clinical, phemap)),
  # [human_review]: check ICD audit results 
  tar_target(write_audit_tsv,
    {
      dir.create("scratch", showWarnings = FALSE, recursive = TRUE)  # [FIX] ensure dir
      out_path <- "scratch/embed_ambiguous_icd.tsv"
      data.table::fwrite(audit_results, file = out_path, sep = '\t')
      return(out_path) # Returns path string so targets knows it ran successfully
    },
    format = "file"
  ),
  # [human_review]: eyeball Excel spreadsheet of ICD audit 
  tar_target(write_audit_xlsx,
    {
      dir.create("data/review", showWarnings = FALSE, recursive = TRUE)  # [FIX] ensure dir
      out_path <- "data/review/embed_ambiguous_icd.xlsx"
      WriteAuditXlsx(audit_results, out_path)
      return(out_path) # Returns path string so targets knows it ran successfully
    },
    format = "file"
  ),
  # [human_review]: inspect QC range-code typing (V-codes -> ambiguous,
  #       ICD-9 procs 00-99.99 -> not_icd_cm). 
  tar_target(qc_range_codetypes,
    {
      tmp <- unique(embed_clinical[is_range == TRUE,
                    .(code_type, range_min, range_max, ambiguous)])
      setorder(tmp, code_type, range_min, range_max, ambiguous)
      dir.create("scratch", showWarnings = FALSE, recursive = TRUE)
      out_path <- "scratch/qc_range_codetypes.tsv"
      fwrite(tmp, out_path, sep = "\t")
      return(out_path)
    },
    format = "file"
  ),
  # add missing parent nodes to KESER
  tar_target(dict_keser_full, ExpandMissingParents(dict_wiki_keep, embed_clinical)),
  # query phecodes to validate 
  #       sort() drops the NA query phecodes (COVID, Long_COVID: no canonical
  #       1.2 code) via na.last = NA, which is the intended behaviour.
  tar_target(phecodes_to_validate, embed_clinical[, sort(unique(phecode_query))]),
  # sort edges so u < v for all edges(u,v)
  # THEN restrict to the query-anchored bipartite projection (U_subset x V_full).
  tar_target(dt_true,
    CanonicalizeEdges(dict_keser_full[, .(u = phe1, v = phe2)])[
      u %chin% phecodes_to_validate | v %chin% phecodes_to_validate]
  ),
  # per-query synthetic controls: deterministic set filtering, then sampling. 
  # BuildValidUniverse() filters out true edges (dt_true) and 
  #   hierarchy pairs (ancestor/descendant) pairs 
  # strata sample from the set of possible edges after filtering
  # except degree_matched, which perturbs the non-query endpoint of true edges. 
  # Details in ./R/functions.R.
  tar_target(dt_synthetic_controls,
    GenerateSyntheticControls(dt_true, phe_nodes, phecodes_to_validate, seed = 20260528)
  ),

  # ---- Evaluation stage ----------------------------------------------------
  # Predicted scores, canonicalized (u>v) then aggregated (max/mean/min per pair).
  tar_target(dt_predicted,
    AggregatePredicted(embed_clinical, u_col = "phecode_query", v_col = "phecode")),
  # Labelled eval set: positives broadcast into every (stratum, grid_size).
  tar_target(dt_eval, BroadcastPositives(dt_true, dt_synthetic_controls)),
  # Attach scores, long over score_type, floor unscored pairs.
  tar_target(dt_scored, ScoreEdges(dt_eval, dt_predicted, unscored = "floor")),
  # Scalar metrics per stratum x grid_size x score_type.
  tar_target(dt_auc, ComputeAUC(dt_scored)),
  # Full threshold sweep (sens/spec/PPV/NPV). LARGE (stratums x grid_sizes)
  tar_target(dt_roc, ComputeOperatingPoints(dt_scored)),
  # Plot object.
  tar_target(plot_auc, PlotAUC(dt_auc)),
  # Per-query plotting metadata
  tar_target(query_phecode_colors,
    c(
      "411.4"  = "#E69F00",  # CAD
      "496.21" = "#56B4E9",  # COPD
      "296.2"  = "#009E73",  # Depression
      "389"    = "#F0E442",  # Hearing_Loss
      "428"    = "#0072B2",  # HF
      "411.2"  = "#D55E00",  # MI
      "335"    = "#CC79A7",  # MS
      "714.1"  = "#999999",  # RA
      "250.2"  = "#000000"   # T2DM
    )
  ),
  tar_target(query_phecode_labels,
    c(
      "411.4"  = "CAD",
      "496.21" = "COPD",
      "296.2"  = "Depression",
      "389"    = "Hearing Loss",
      "428"    = "HF",
      "411.2"  = "MI",
      "335"    = "MS",
      "714.1"  = "RA",
      "250.2"  = "T2DM"
    )
  ),
  # ---- Outputs (tracked as files) -----------------------------------------
  # Core sinks (auc table, roc table, auc-vs-gridsize plot) go through the
  # single WriteEvalOutputs() writer, which creates OUT_DIR and returns the
  # three file paths so targets tracks them as a set.
  tar_target(out_core,
    WriteEvalOutputs(dt_auc, dt_roc, plot_auc, outdir = OUT_DIR),
    format = "file"),
  # Per-query AUC bar charts (random stratum, score_mean, 1:1). Both AUROC and
  # AUPRC are written, per the project rule to always report them together.
  tar_target(out_auc_by_query,
    {
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
      f_roc <- file.path(OUT_DIR, "auc_by_query_auroc.png")
      f_prc <- file.path(OUT_DIR, "auc_by_query_auprc.png")
      ggsave(f_roc, PlotAUCByQuery(dt_auc, metric = "auroc",
             palette = query_phecode_colors, labels = query_phecode_labels),
             width = 11, height = 8, dpi = 300)
      ggsave(f_prc, PlotAUCByQuery(dt_auc, metric = "auprc",
             palette = query_phecode_colors, labels = query_phecode_labels),
             width = 11, height = 8, dpi = 300)
      c(f_roc, f_prc)
    },
    format = "file"
  ),
  # Per-query operating-characteristic panels: one PNG per query phecode
  # (2x2 sensitivity/specificity/PPV/NPV vs threshold). Returns the vector of
  # written paths so targets tracks every file.
  tar_target(out_roc_by_query,
    {
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
      plots <- PlotOperatingPointsByQuery(dt_roc, labels = query_phecode_labels)
      paths <- vapply(names(plots), function(q) {
        f <- file.path(OUT_DIR, sprintf("roc_query_%s.png", q))
        ggsave(f, plots[[q]], width = 8, height = 7, dpi = 300)
        f
      }, character(1))
      unname(paths)
    },
    format = "file"
  ),
  # Headline per-query AUC summary table (random stratum, score_mean, 1:1).
  tar_target(out_headline_auc,
    {
      dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
      f <- file.path(OUT_DIR, "headline_auc_random_1to1.csv")
      fwrite(HeadlineAucTable(dt_auc, labels = query_phecode_labels), f)
      f
    },
    format = "file"
  )
)
