# Custom QC functions ####
## LOOKS >= ${threshold} CHARACTERS (char) IN COLUMN (col_name) of data.table (dt)
# Check for >= 2 hyphens:   qc_pattern_count(code_ranges, "code", "-", 2)
# Check for >=3 semicolons: qc_pattern_count(code_ranges, "code", ";", 3)
QcPatternCount <- function(dt, col_name, char, threshold = 2) {
  # Logic Summary:
  # 1. Constructs regex to find 'n' occurrences of 'char' using quantifier {n,}.
  # 2. Uses fixed = FALSE to allow regex metacharacter escaping if needed.
  # 3. Halts execution via stop() if nrow > 0.
  pattern <- paste0("(", char, ".*){", threshold, ",}")
  tryCatch({
    problem_rows <- dt[grepl(pattern, get(col_name))]
    row_count <- nrow(problem_rows)
    
    if (row_count > 0) {
      stop(sprintf("%d rows with >= %d occurrences of '%s'", 
                   row_count, threshold, char))
    }
    return(sprintf("QC PASS: %d rows with >= %d occurrences of '%s'", row_count, threshold, char))
  }, error = function(e) {
    stop(paste0("QC FAIL: ", e$message))
  })
}


LoadPhecodeMap <- function(path) {
  phemap <- fread(path, colClasses = "character")
  
  cols_all  <- c("ICD", "Phecode", "ICDString", "PhecodeString", "PhecodeCategory")
  cols_base <- c("ICD", "Phecode")
  
  row_diff <- uniqueN(phemap, by = cols_all) - uniqueN(phemap, by = cols_base)
  if (row_diff > 0) {
    warning("Phecode-ICD mapping table has duplicated rows: ", row_diff, " extra rows")
  }
  
  phemap <- unique(phemap[, ..cols_all])
  phemap[, root := tstrsplit(ICD, ".", fixed = TRUE)[[1]]]
  
  # Explicit return of materialized data.table (avoids returning lazy.frame)
  return(phemap)
}

# Data Loading Functions ####
PopulateDictionary <- function(file_path, type = "ICD", subtype = "diagnosis", source = "OHDSI") {
  if (type == "ICD" && subtype == "diagnosis" && source == "OHDSI") {
    dict <- data.table::fread(
      file = file_path, 
      sep = "\t", 
      quote = "", 
      select = c("concept_code", "vocabulary_id"))[
        # select only ICD-CM codes, no Procs or WHO
        vocabulary_id %chin% c("ICD9CM", "ICD10CM")
      ]
    # de-dup in case input table wasn't
    dict <- unique(dict, by = c("concept_code", "vocabulary_id"))
    # flag codes used in multiple voacbularies
    # e.g. V-codes in ICD9 & ICD10)
    ambiguous_codes <- dict[, .N, by = concept_code][N > 1, concept_code] 
    # add ambiguous column to dictionary for downstream processing
    dict[, ambiguous := concept_code %chin% ambiguous_codes]
    # output cols
    return(dict)
  }
  stop("Dictionary type/source combination not implemented.")
}

LoadKeserDxDx <- function(filepath) {
  temp_env <- new.env(parent = emptyenv())  # emptyenv() is lighter — no need to inherit base
  load(
    filepath,                                # [FIX] was a hardcoded path; honour the arg
    envir = temp_env
  )
  dt <- as.data.table(temp_env$dict.wiki.keep)[
                      phecode1 != phecode2
                    ][
                      , src := "keser"
                    ]
  # remove self-references
  rm(temp_env)
  # 1a) strip KESER prefix "PheCode:"
  # input:
  #phecode1	phecode2	src
  #PheCode:540.1	PheCode:442.1	keser
  #PheCode:594.1	PheCode:442.1	keser
  # output:
  # phe1 phe2
  # 540.1 442.1
  dt[, c("phe1", "phe2") := lapply(.SD, function(x) 
          stringi::stri_replace_first_fixed(x, "PheCode:", "")),
          .SDcols = c("phecode1", "phecode2")
        ]
  return(dt)
}

LoadEmbedDx <- function(directory) {
  fnames <- list.files(path = directory, pattern = "_Diagnosis\\.csv$")
  
  flist <- lapply(file.path(directory, fnames), fread)
  names(flist) <- gsub("_Diagnosis\\.csv$", "", fnames)
  dt = rbindlist(flist, idcol = "clinical_query")

  # map clinical query to phecode 
  phecode_1_2_map <- c(
    "CAD"          = "411.4",
    "COPD"         = "496.21",
    "COVID"        = NA_character_, # Not in canonical 1.2
    "Depression"   = "296.2",
    "Hearing_Loss" = "389",
    "HF"           = "428",
    "Long_COVID"   = NA_character_, # Not in canonical 1.2
    "MI"           = "411.2",
    "MS"           = "335",
    "RA"           = "714.1",
    "T2DM"         = "250.2"
  )
  dt[, phecode_query := unname(phecode_1_2_map[clinical_query])]
  return(dt)
}

# Data Munging functions ####
# Pre-processes data.table (DT) with ICD code column (col_name)
SplitRanges <- function(DT, col_name = "code") {
  # remove garbage rows; need individual codes or code ranges, not entire ICD dictionary
  DT <- DT[target != 'ICD-10-CM', env = list(target = col_name)]
  # embed_clinical format:
  #   rank,code,description,score_cosine,clinical_query
  #   1,U07.1,covid-19,1.0,COVID
  #   1,714.0,rheumatoid arthritis,1.0,RA
  #   1,"340,G35",multiple sclerosis,0.99,MS
  # Split rows with compound codes (e.g. "411,412") into multiple rows
  # env = evaluate col_name from the function input only
  DT <- DT[, stringi::stri_split_fixed(target, ","), 
            by = setdiff(names(DT), col_name),
            env = list(target = col_name)]
  setnames(DT, 'V1', col_name)
  # Add flag for whether code is a code range
  # if code is single code, range_min = code
  range_cols = c('range_min', 'range_max')
  DT[, (range_cols) := tstrsplit(target, "-", fixed = TRUE, keep = 1:2), 
      env = list(target = col_name)]
  # ICDs only use uppercase letters and strip whitespace 
  code_cols = c(col_name, 'range_min', 'range_max')
  DT[, (code_cols) := lapply(.SD, function(x) toupper(trimws(x))), .SDcols = code_cols]
  DT[, is_range := !is.na(range_max)]
  # split the ranges into root code and digits after decimal
  DT[, c("min_root", "min_decimal") := tstrsplit(range_min, ".", fixed = TRUE)]
  DT[is_range == TRUE, 
    c("max_root", "max_decimal") := tstrsplit(range_max, ".", fixed = TRUE)]  
  return(DT)
}

# non-equi joins for character vectors in R data.table
JoinStringRanges <- function(dict_dt, query_dt, dict_col, min_col, max_col, return_cols, arg_mult = "first") {
  
  # 1. Unified Radix Sorting
  all_strings <- sort(
    unique(c(dict_dt[[dict_col]], query_dt[[min_col]], query_dt[[max_col]])),
    method = "radix"
  )
  
  # 2. Isolate dictionary and calculate taxonomy dimensions
  target_cols <- c(dict_col, return_cols)
  d_temp <- dict_dt[, .SD, .SDcols = target_cols]
  
  d_temp[, `:=`(pt_int  = match(TARGET, all_strings),
                str_len = nchar(TARGET)),
         env = list(TARGET = dict_col)]
  
  setorder(d_temp, pt_int)
  
  # 3. Map queries into integer space WITH tracking ID
  q_temp <- data.table(
    q_id    = seq_len(nrow(query_dt)), # Unique identifier for expansion tracking
    min_int = match(query_dt[[min_col]], all_strings),
    max_int = match(query_dt[[max_col]], all_strings),
    q_len   = nchar(query_dt[[min_col]]) 
  )
  
  # Prefix return columns
  ret_names <- paste0("ret_", return_cols)
  setnames(d_temp, return_cols, ret_names)
  
  # 4. Multi-condition Binary Merge
  matched <- d_temp[
    q_temp, 
    on = .(str_len == q_len, pt_int >= min_int, pt_int <= max_int), 
    mult = arg_mult,
    allow.cartesian = TRUE # Required for mult="all" interval expansions
  ]
  
  # Return the tracking ID alongside the payload
  return(matched[, .SD, .SDcols = c("q_id", ret_names)])
}

# Check if ICD code is ICD-9, ICD-10, or ambiguous
GetCodeType <- function(DT, col_to_check, dict_dt, type = "ICD") {
  
  # initialize columns 
  if (!"code_type" %in% names(DT)) { DT[, code_type := NA_character_] }
  if (!"ambiguous" %in% names(DT)) { DT[, ambiguous := NA] }
  
  # check for ranges
  if (any(stringi::stri_detect_regex(DT[[col_to_check]], "[,\\-]"))) {
    stop(sprintf("Detected unsplit ranges/lists in `%s`. Run SplitRanges() first.",
    col_to_check))
  }
  
  # checks dictionary loaded by PopulateDictionary()
  if (type == "ICD") {
    # Primary join: exact match
    DT[dict_dt, on = setNames("concept_code", col_to_check),
       `:=`(code_type = i.vocabulary_id, 
            ambiguous = i.ambiguous)]

    # Fallback: base codes (e.g. E849) absent from dict but whose sub-codes
    # (e.g. E849.1) are present. Derive a root-level lookup from dict.
    dict_roots <- dict_dt[
      # only rows that actually have a decimal (i.e. have a root)
      grepl(".", concept_code, fixed = TRUE),
      .(root = tstrsplit(concept_code, ".", fixed = TRUE)[[1]], 
        vocabulary_id, ambiguous)
    ][,
      .(vocabulary_id = if (uniqueN(vocabulary_id) == 1) vocabulary_id[1] 
                        else "Ambiguous",
        ambiguous     = any(ambiguous) | uniqueN(vocabulary_id) > 1
      ),
      by = root
    ]

    # merging root codes that dictionary didn't define explicitly
    DT[dict_roots, on = setNames("root", col_to_check),
       `:=`(code_type = fcoalesce(code_type, i.vocabulary_id), 
            ambiguous = fcoalesce(ambiguous, i.ambiguous))]
    DT[is.na(code_type), `:=`(code_type = "not_icd_cm", ambiguous = FALSE)]

    # 1. Isolate the ranges that failed both primary and root point-lookups.
    # 2. Execute an O(log n) non-equi join against dict_roots using the pre-parsed min/max boundaries.
    # 2b. data.table can't do character non-equi joins
    # 2c. create a radix dictionary
    range_mask <- DT[, is_range == TRUE & code_type == "not_icd_cm"]

    if (any(range_mask, na.rm = TRUE)) {
      # generalized string range join
      matched <- JoinStringRanges(
        dict_dt     = dict_roots, 
        query_dt    = DT[(range_mask)], 
        dict_col    = "root", 
        min_col     = "min_root", 
        max_col     = "max_root", 
        return_cols = c("vocabulary_id", "ambiguous")
      )
      
      # update code_type for code ranges that weren't classified
      DT[(range_mask), `:=`(
        code_type = fcoalesce(matched$ret_vocabulary_id, "not_icd_cm"),
        ambiguous = fcoalesce(matched$ret_ambiguous, FALSE)
      )]
    }

    return(DT)
  } # end ICD code typing
  stop("Validation type/ontology not implemented.")
}

MapIcdToPhecode <- function(icd_dt, phemap) {
  # Deep copy prevents in-place mutation of the upstream targets cache
  dt <- data.table::copy(icd_dt)
  dt[, row_id := .I]
  
  # PART 1: EXACT MATCH (Point Codes)
  matched_exact <- phemap[
    dt[!ambiguous & !is_range], 
    on = .(ICD = code), 
    nomatch = NULL,
    .(row_id, Phecode, ICDString, PhecodeString, PhecodeCategory)
  ]
  
  dt[matched_exact, on = .(row_id), `:=`(
    phecode          = i.Phecode,
    icd_string       = i.ICDString,
    phecode_string   = i.PhecodeString,
    phecode_category = i.PhecodeCategory
  )]
  
  # Clean up global row tracking
  dt[, row_id := NULL]
  
  # ==========================================
  # PART 2: INTERVAL EXPANSION (Range Codes)
  # ==========================================
  range_mask <- dt[, !ambiguous & is_range]
  
  if (any(range_mask, na.rm = TRUE)) {
    # 1. Isolate the subset and assign the tracking ID to match q_temp
    ranges_to_expand <- dt[(range_mask)][, q_id := .I]
    # 2. Inject "ICD" into the return_cols array
    matched_ranges <- JoinStringRanges(
      dict_dt     = phemap,
      query_dt    = ranges_to_expand,
      dict_col    = "root",
      min_col     = "min_root",
      max_col     = "max_root",
      return_cols = c("ICD", "Phecode", "ICDString", "PhecodeString", "PhecodeCategory"), 
      arg_mult    = "all"
    )
    # 2. Map ret_ICD to a dedicated QC column to preserve the original `code` column
    expanded_dt <- matched_ranges[ranges_to_expand, on = "q_id", nomatch = NULL]
    expanded_dt[, `:=`(
      matched_icd      = ret_ICD,
      phecode          = ret_Phecode,
      icd_string       = ret_ICDString,
      phecode_string   = ret_PhecodeString,
      phecode_category = ret_PhecodeCategory
    )]
    
    # 3. Vectorized Sieve
    expanded_dt <- expanded_dt[
      # Appending "Z" : "I29.9" <= "I29Z" = True
      matched_icd >= range_min & matched_icd <= paste0(range_max, "Z")
    ]
    # 4. remove temp cols
    temp_cols <- c("q_id", paste0("ret_", 
      c("ICD", "Phecode", "ICDString", "PhecodeString", "PhecodeCategory")))
    expanded_dt[, (temp_cols) := NULL]
    
    # 5. Reconstruct Master Table
    dt <- rbindlist(list(
      dt[!(range_mask)], # mapped single codes
      expanded_dt       # mapped code ranges
    ), use.names = TRUE, fill = TRUE)
  }
  return(dt)
}

# Try to map ambiguous ICD codes (V-block) to Phecode
GenerateAmbiguityAudit <- function(clinical_dt, phemap) {
  # 1. Isolate Ambiguous Point Codes
  audit_dt <- clinical_dt[ambiguous & !is_range, .(code, description)]
  
  # 2. String Normalization
  # Collapses multiple spaces and drops all punctuation to maximize exact match yield
  NormalizeString <- function(v) {
    tolower(trimws(gsub("[[:punct:]]", "", gsub("\\s+", " ", v))))
  }
  
  audit_dt[, norm_desc := NormalizeString(description)]
  # [FIX] copy before adding a column; otherwise norm_icd_str is written into
  #       the caller's phemap (a shared targets object) by reference
  phemap <- copy(phemap)
  phemap[, norm_icd_str := NormalizeString(ICDString)]
  
  # Exact String Join
  audit_matched <- phemap[
    audit_dt, 
    on = .(ICD = code, norm_icd_str = norm_desc),
    .(code        = i.code, 
      embed_desc = i.description, 
      Phecode     = x.Phecode, 
      dict_desc   = x.ICDString)
      #, nomatch = NULL
  ]
  
  # 3. Match Statistics
  match_rate <- sum(!is.na(audit_matched$Phecode)) / nrow(audit_matched)
  message(sprintf("Exact String Match Rate: %.2f%%", match_rate * 100))
  
  # 4. Process Unmatched Codes
  unmatched <- audit_matched[is.na(Phecode), .(code, embed_desc)]
  
  # Generate Cartesian expansion mapping the unmatched code to ALL possible dictionary definitions
  eval_dt <- phemap[
    unmatched, 
    on = .(ICD = code), 
    allow.cartesian = TRUE,
    .(code = i.code, embed_desc = i.embed_desc, 
      dict_desc = x.ICDString, Phecode = x.Phecode)
  ]

  # Isolate orphan codes entirely missing from the dictionary
  orphan_mask <- is.na(eval_dt$Phecode)
  orphans <- eval_dt[(orphan_mask), .(code, embed_desc)]
  
  # Drop orphans from the string evaluation to prevent NA propagation
  eval_dt <- eval_dt[!(orphan_mask)]
  
  # 5. String Similarity Evaluation
  # Cosine similarity of character 4-grams
  eval_dt[, sim_score := stringdist::stringsim(
    embed_desc, dict_desc, method = "cosine", q = 4
  )]
  # [FIX] previously took .SD[1] BEFORE sorting, so it kept an arbitrary dict
  #       match per code rather than the closest one. Order within code by
  #       descending similarity (NA sims pushed last) THEN keep the top row.
  setorder(eval_dt, code, -sim_score, na.last = TRUE)
  eval_dt = eval_dt[, .SD[1], by = code]
  # present highest-similarity matches first (original behaviour); flip the
  # sign if you'd rather surface the least-confident codes for review triage
  setorder(eval_dt, -sim_score, na.last = TRUE)

  eval_dt[, `:=`(
  mismatch     = NA_character_, # Blank column for mismatch
  human_review = NA_character_  # Blank column for clinical justification
)]
  
  # Extract the single closest dictionary match for each unresolved code
  return(eval_dt)
}

ProcessClinicalData <- function(predict_dir, dict_icd, phemap) {
  dt <- LoadEmbedDx(predict_dir)
  # split lists (code1,code2,etc) and ranges (eg R10-R19) into multiple rows
  dt <- SplitRanges(dt, "code")
  # classify code as ICD9CM, ICD10CM, ambiguous, NA
  dt <- GetCodeType(dt, "range_min", dict_icd, type = "ICD")
  # Map ICD codes to Phecode using Phecode Mapping 1.2 (NO ICD-10 U-codes)
  dt <- MapIcdToPhecode(dt, phemap)
  return(dt)
}

WriteAuditXlsx <- function(audit_results, out_path) {
      wb <- openxlsx::createWorkbook()
      # initialize worksheet
      openxlsx::addWorksheet(wb, "Ambiguous_ICDs")
      # Write the data
      openxlsx::writeData(wb, sheet = "Ambiguous_ICDs", x = audit_results, 
        headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
      )
      # Freeze the top row and adjust column widths for human readability
      openxlsx::freezePane(wb, sheet = "Ambiguous_ICDs", firstRow = TRUE)
      openxlsx::setColWidths(wb, sheet = "Ambiguous_ICDs", 
                       cols = 1:6, widths = c(10, 50, 50, 12, 20, 30))
      # Save to the review directory
      openxlsx::saveWorkbook(wb, file = out_path, overwrite = TRUE)
}


# Hierarchical knowledge graph expansion ####
## Example: check_phecode_hierarchy(DT = dict_wiki_keep)
CheckPhecodeHierarchy <- function(DT) {
  edges <- DT[, .(
    id = .I,
    p1 = sub("PheCode:", "", phecode1),
    p2 = sub("PheCode:", "", phecode2)
  )]

  get_lineage <- function(dt, col_name) {
    lineage <- dt[, .(
      id   = id,
      leaf = get(col_name),
      root = sub("\\..*", "", get(col_name)),
      mid  = fifelse(grepl("\\.[0-9]{2}$", get(col_name)), 
        sub(".$", "", get(col_name)), NA_character_)
    )]
    melt(lineage, id.vars = "id", measure.vars = c("leaf", "root", "mid"),
         value.name = "node", na.rm = TRUE)[, .(id, node)]
  }

  p1_lineage <- get_lineage(edges, "p1")
  p2_lineage <- get_lineage(edges, "p2")

  theoretical_edges <- unique(
    p1_lineage[p2_lineage, on = "id", allow.cartesian = TRUE][, .(p1 = node, p2 = i.node)]
  )
  empirical_edges <- unique(edges[, .(p1, p2)])
  missing_edges   <- theoretical_edges[!empirical_edges, on = .(p1, p2)]

  n_theoretical <- nrow(theoretical_edges)
  n_missing     <- nrow(missing_edges)
  compliance    <- round(100 * (n_theoretical - n_missing) / n_theoretical, 2)

  if (n_missing > 0) {
    warning(sprintf(
      "Phecode hierarchy not fully closed: %d missing edges (compliance %.2f%%)",
      n_missing, compliance
    ))
  }

  return(data.table(
    metric = c("empirical_edges", "theoretical_edges", "missing_edges", "compliance_pct"),
    value  = c(nrow(empirical_edges), n_theoretical, n_missing, compliance)
  ))
}


#' Build long-form lineage table: each node -> each of its ancestor codes.
#'
#' Uses character truncation of the post-decimal segment to generate the full
#' ancestor chain. For "S52.521A" this yields c("S52", "S52.5", "S52.52",
#' "S52.521"). Works uniformly for PheCodes, ICD-9, and ICD-10-CM (including
#' 7th-character extensions). Codes with no decimal point have no ancestors.
#'
#' Note: assumes at most one decimal per code (true for PheCodes/ICD-9/ICD-10-CM).
#' Malformed inputs with multiple dots will put extra dots into the "post"
#' segment and produce nonsense ancestors that simply won't match anything.
build_lineage <- function(nodes) {
  dt <- data.table(node = unique(nodes))
  dt[, has_dot := grepl(".", node, fixed = TRUE)]
  dt <- dt[has_dot == TRUE]
  if (nrow(dt) == 0L) {
    return(data.table(node = character(0), parent = character(0)))
  }

  dt[, pre  := sub("\\..*", "", node)]
  dt[, post := sub("^[^.]*\\.", "", node)]
  dt <- dt[nchar(post) > 0L]
  if (nrow(dt) == 0L) {
    return(data.table(node = character(0), parent = character(0)))
  }
  dt[, post_len := nchar(post)]

  expanded <- dt[, .(level = seq_len(post_len)), by = .(node, pre, post)]
  expanded[, parent := fifelse(
    level == 1L,
    pre,
    paste0(pre, ".", substring(post, 1L, level - 1L))
  )]

  expanded[node != parent, .(node, parent)]
}

#' Expand a KESER co-occurrence graph to cover missing retrieved parent phecodes.
#'
#' The embedding model embeds a clinical query (e.g. "Rheumatoid Arthritis", "CAD") and
#' returns matched phecodes via `embed_clinical`. Some of those phecodes are
#' higher-level parents that are absent from the KESER co-occurrence graph
#' (`dict_wiki_keep`), even though their child phecodes are present. For each
#' empirical KESER edge (A, B) whose endpoint has a missing-parent ancestor
#' A', this function generates the inherited edge (A', B), letting downstream
#' graph queries reach the missing parents via their children's neighborhoods.
#'
#' The set of missing parents is computed internally as
#'   setdiff(embed_clinical[ambiguous == FALSE, unique(phecode)],
#'           unique(c(dict_wiki_keep$phe1, dict_wiki_keep$phe2)))
#' then filtered to only those that are actually ancestors of some KESER node
#' (per `build_lineage`). The `ambiguous == FALSE` filter is hardcoded — only
#' unambiguous retrieved phecodes contribute.
#'
#' Lineage is parsed via character truncation of the post-decimal segment, so
#' the function handles arbitrary phecode depth (PheCodes, and by extension
#' ICD-9 / ICD-10-CM if the input ever broadens).
#'
#' @param dict_wiki_keep data.table with columns `phe1`, `phe2` (plus any
#'        optional metadata like cosine/weight). The KESER edge list.
#' @param embed_clinical data.table with columns `phecode` and `ambiguous`
#'        (logical). Only rows with `ambiguous == FALSE` are used.
#' @param undirected if TRUE (default), treat (A, B) and (B, A) as the same
#'        edge. Canonicalizes both theoretical and empirical edges to
#'        (min, max) before the anti-join. Set FALSE if the graph is directed.
#' @param keep_self_loops if FALSE (default), drop expanded edges with
#'        phe1 == phe2. These can arise when both endpoints of an empirical
#'        edge share the same missing ancestor.
#' @return `dict_wiki_keep` with new rows appended (`expanded = TRUE`) and an
#'         `expanded` flag on all rows (`FALSE` on original rows). Non-edge
#'         metadata on new rows is NA — downstream weighting must handle this.
ExpandMissingParents <- function(dict_wiki_keep, 
                                  embed_clinical, 
                                  undirected = TRUE, 
                                  keep_self_loops = FALSE) {
  stopifnot(
    is.data.table(dict_wiki_keep),
    all(c("phe1", "phe2") %in% names(dict_wiki_keep)),
    is.data.table(embed_clinical),
    all(c("phecode", "ambiguous") %in% names(embed_clinical))
  )

  # ---- 0. Identify missing parent codes -------------------------------------
  embed_phecodes <- embed_clinical[ambiguous == FALSE, 
                  unique(c(phecode_query, phecode))]
  keser_phecodes  <- unique(c(dict_wiki_keep$phe1, dict_wiki_keep$phe2))
  missing_codes   <- setdiff(embed_phecodes, keser_phecodes)

  if (length(missing_codes) == 0L) {
    message("All unambiguous retrieved phecodes are already in KESER; no expansion needed.")
    out <- copy(dict_wiki_keep)
    if (!"expanded" %in% names(out)) out[, expanded := FALSE]
    return(out[])
  }

  # ---- 1. Lineage generation ------------------------------------------------
  edges        <- dict_wiki_keep[, .(id = .I, phe1, phe2)]
  unique_nodes <- unique(c(edges$phe1, edges$phe2))
  lineage_long <- build_lineage(unique_nodes)

  # ---- 2. Intersection test -------------------------------------------------
  # Keep only missing codes that are actually ancestors of some KESER node.
  # A missing retrieved code with no children in KESER cannot be reached by
  # neighborhood inheritance — flag it but don't try to expand.
  valid_missing_parents   <- intersect(missing_codes, lineage_long$parent)
  unreachable_missing     <- setdiff(missing_codes, valid_missing_parents)

  if (length(unreachable_missing) > 0L) {
    message(sprintf(
      "%d missing retrieved phecode(s) have no children in KESER and cannot be expanded: %s",
      length(unreachable_missing),
      paste(unreachable_missing, collapse = ", ")
    ))
  }

  if (length(valid_missing_parents) == 0L) {
    message("No missing retrieved phecodes are ancestors of any KESER node; nothing to expand.")
    out <- copy(dict_wiki_keep)
    if (!"expanded" %in% names(out)) out[, expanded := FALSE]
    return(out[])
  }

  # ---- 3. Expansion map (node -> {node itself, qualifying ancestors}) -------
  node_map   <- data.table(node = unique_nodes, exp_node = unique_nodes)
  parent_map <- lineage_long[parent %in% valid_missing_parents,
                             .(node, exp_node = parent)]
  full_map   <- rbind(node_map, parent_map)

  # ---- 4. Cartesian edge generation by original edge id ---------------------
  exp_p1 <- full_map[edges, on = .(node = phe1),
                     .(id, p1_exp = exp_node), allow.cartesian = TRUE]
  exp_p2 <- full_map[edges, on = .(node = phe2),
                     .(id, p2_exp = exp_node), allow.cartesian = TRUE]

  theoretical_edges <- unique(
    exp_p1[exp_p2, on = "id", allow.cartesian = TRUE][,
      .(phe1 = p1_exp, phe2 = p2_exp)]
  )

  # ---- 5. Canonicalize if undirected, then anti-join ------------------------
  empirical_edges <- unique(edges[, .(phe1, phe2)])

  if (undirected) {
    canonicalize <- function(d) {
      d[, `:=`(phe1 = pmin(phe1, phe2),
               phe2 = pmax(phe1, phe2))][]
    }
    canonicalize(theoretical_edges)
    canonicalize(empirical_edges)
    theoretical_edges <- unique(theoretical_edges)
    empirical_edges   <- unique(empirical_edges)
  }

  new_edges <- theoretical_edges[!empirical_edges, on = .(phe1, phe2)]

  if (!keep_self_loops) {
    new_edges <- new_edges[phe1 != phe2]
  }

  # ---- 6. Flag and bind back ------------------------------------------------
  new_edges[, expanded := TRUE]
  out <- copy(dict_wiki_keep)
  if (!"expanded" %in% names(out)) out[, expanded := FALSE]

  result <- rbind(out, new_edges, fill = TRUE)

  message(sprintf(
    "Expanded KESER graph by %d edges via %d missing retrieved parent phecode(s).",
    nrow(new_edges), length(valid_missing_parents)
  ))

  return(result[])
}

# -------------------------------------------------------------------------
# Helper Functions for Graphs
# -------------------------------------------------------------------------
# Edge Canonicalization Function
CanonicalizeEdges <- function(dt, col1 = "u", col2 = "v") {
  dt <- copy(dt)  # [FIX] avoid mutating the caller's data.table by reference
  # Use distinct macro tokens (C1, C2) to prevent AST symbol collisions
  dt[, c("u", "v") := .(
      fifelse(C1 < C2, C1, C2),
      fifelse(C1 < C2, C2, C1)
    ), env = list(C1 = col1, C2 = col2)
  ] 
  # Filter out trivial self-edges
  return(unique(dt[u != v]))
}

# check if Phecode node u is an ancestor of Phecode node v
IsAncestorDescendant <- function(u, v) {
  # u: Character vector (lexicographically first, potential ancestor)
  # v: Character vector (lexicographically second, potential descendant)
  # PheCode mechanics: 'v' is a descendant if it strictly extends the string 'u'
  return(startsWith(v, u) & (nchar(v) > nchar(u)))
}

# -------------------------------------------------------------------------
# Candidate Universe (built ONCE, shared across enumerable strata)
# -------------------------------------------------------------------------
# Enumerates every canonical query-anchored pair {u, v} with u in the query
# set, strips positives (dt_true) and ancestor/descendant pairs, and attaches
# category numbers. The phecode graph is small enough (~1.8k nodes) that the
# query x full candidate space is trivially enumerable, so the random /
# non-adjacent / intracategory strata become filters + sampling WITHOUT
# replacement rather than rejection loops.
#
# Assumes dt_true is canonicalized (u < v) -- enforced in GenerateSyntheticControls.
BuildValidUniverse <- function(query_phecodes, dt_nodes, dt_true) {
  all_phecodes <- dt_nodes$phecode

  # CJ cross-join; query-query pairs appear as (a,b) and (b,a) and collapse
  # to a single canonical row after fifelse + unique.
  univ <- CJ(a = query_phecodes, b = all_phecodes, unique = TRUE)[a != b]
  univ[, `:=`(u = fifelse(a < b, a, b),
              v = fifelse(a < b, b, a))]
  univ <- unique(univ[, .(u, v)])

  univ <- univ[!dt_true, on = .(u, v)]          # drop true edges (positives)
  univ <- univ[!IsAncestorDescendant(u, v)]     # drop hierarchy pairs

  # Categories attached once; non-adjacent / intracategory just filter on these.
  univ[dt_nodes, cat_u := i.category_number, on = .(u = phecode)]
  univ[dt_nodes, cat_v := i.category_number, on = .(v = phecode)]
  return(univ)
}

# -------------------------------------------------------------------------
# Per-query sampling helpers
# -------------------------------------------------------------------------
# Negative controls are matched to positives PER QUERY PHECODE, so grid_size is
# the negative:positive RATIO (1, 2, 5, ...) rather than an absolute count. A
# global absolute count cannot be 1:1 for every phenotype at once because each
# phenotype has a different number of positive edges; anchoring per query makes
# the realized ratio equal grid_size (until a query's candidate pool runs out).

# Explode an edge/candidate table to one row per QUERY endpoint. A pair whose u
# and v are BOTH query phecodes appears twice (once per query), so it is counted
# under each query's evaluation. Other columns (e.g. cat_u/cat_v) carry through.
AnchorByQueryEndpoint <- function(dt, query_phecodes) {
  on_u <- dt[u %chin% query_phecodes]
  on_u[, phecode_query := u]
  on_v <- dt[v %chin% query_phecodes]
  on_v[, phecode_query := v]
  rbindlist(list(on_u, on_v), use.names = TRUE)
}

# Slice an already-shuffled (u, v) pool into nested ratio cells: head(pool,
# n_pos * r) over increasing r yields nested subsets (smaller ratios subset
# larger ones). Each returned cell carries phecode_query, stratum, grid_size
# (= ratio) and true_edge = 0.
SliceRatioCells <- function(ordered, n_pos, ratios, stratum_label, q) {
  n_have <- nrow(ordered)
  rbindlist(lapply(ratios, function(r) {
    cell <- head(ordered, min(n_pos * r, n_have))[, .(u, v)]
    cell[, `:=`(phecode_query = q, stratum = stratum_label,
                grid_size = as.integer(r), true_edge = 0L)]
    cell
  }), use.names = TRUE)
}

# Enumerable strata (random / non_adjacent / intracategory): shuffle the
# query-anchored candidate pool once, then ratio-slice. Warns (and caps) if the
# pool is too small to reach the largest requested ratio for this query.
SampleQueryRatioCells <- function(pool, n_pos, ratios, stratum_label, q) {
  r_max  <- max(ratios)
  target <- n_pos * r_max
  n_take <- min(nrow(pool), target)
  if (n_take < target)
    warning(sprintf("[%s / query %s] %d candidates < %d (%d:1); higher ratios capped.",
                    stratum_label, q, nrow(pool), target, r_max))
  ordered <- if (n_take > 0L) pool[sample(.N, n_take), .(u, v)] else pool[0L, .(u, v)]
  SliceRatioCells(ordered, n_pos, ratios, stratum_label, q)
}

# Degree-matched for ONE query: resample true edges incident to q, keep the q
# endpoint, perturb the other to a random node. Preserves q's endpoint-incidence
# (degree) by construction. Bounded iteration with early exit on saturation,
# then ratio-sliced like the other strata.
SampleDegreeMatchedQuery <- function(q, seeds_q, all_phecodes, dt_true,
                                     n_pos, ratios, max_iter = 50L) {
  empty <- data.table(u = character(), v = character(),
                      phecode_query = character(), stratum = character(),
                      grid_size = integer(), true_edge = integer())
  target <- n_pos * max(ratios)
  if (nrow(seeds_q) == 0L || target == 0L) return(empty)

  res  <- data.table(u = character(), v = character())
  iter <- 0L
  while (nrow(res) < target && iter < max_iter) {
    iter  <- iter + 1L
    batch <- ceiling((target - nrow(res)) * 1.5)
    d     <- seeds_q[sample.int(nrow(seeds_q), batch, replace = TRUE), .(u, v)]
    # q sits on u -> perturb v; otherwise q sits on v -> perturb u.
    keep_u   <- d$u == q
    new_node <- sample(all_phecodes, nrow(d), replace = TRUE)
    d[, `:=`(u = fifelse(keep_u, u, new_node),
             v = fifelse(keep_u, new_node, v))]
    # Canonicalize, then drop self-loops / positives / hierarchy pairs.
    d[, `:=`(u = fifelse(u < v, u, v), v = fifelse(u < v, v, u))]
    d <- d[u != v]
    d <- d[!dt_true, on = .(u, v)]
    d <- d[!IsAncestorDescendant(u, v)]
    before <- nrow(res)
    res <- unique(rbindlist(list(res, d[, .(u, v)])))
    if (nrow(res) == before) break             # saturated -> stop spinning
  }
  if (nrow(res) < target)
    warning(sprintf("[degree_matched / query %s] %d of %d after %d iter; ratios capped.",
                    q, nrow(res), target, iter))
  SliceRatioCells(head(res, target), n_pos, ratios, "degree_matched", q)
}

# -------------------------------------------------------------------------
# Main Execution Function
# -------------------------------------------------------------------------
# Returns negative controls with columns (u, v, phecode_query, stratum,
# grid_size, true_edge), where grid_size is the per-query negative:positive
# ratio. Deterministic set filtering (positives + hierarchy removed up front by
# BuildValidUniverse) followed by per-query random sampling -- NOT rejection
# sampling. Queries with zero positives are excluded (nothing to match).
GenerateSyntheticControls <- function(dt_true, dt_nodes, query_phecodes, seed = 42,
                                      ratios = c(1L, 2L, 5L, 10L, 25L, 50L, 100L)) {
  set.seed(seed)
  ratios <- sort(unique(as.integer(ratios)))

  # Defensive canonicalization of positives: the anti-joins require u < v.
  dt_true <- copy(dt_true)
  dt_true[, `:=`(u = fifelse(u < v, u, v), v = fifelse(u < v, v, u))]
  dt_true <- unique(dt_true[u != v, .(u, v)])

  all_phecodes <- dt_nodes$phecode

  # Per-query positive counts (query-query positives counted under both ends).
  pos_anchored <- AnchorByQueryEndpoint(dt_true[, .(u, v)], query_phecodes)
  n_pos_tab    <- pos_anchored[, .(n_pos = .N), by = phecode_query]

  no_pos <- setdiff(query_phecodes, n_pos_tab$phecode_query)
  if (length(no_pos))
    message(sprintf("No positive edges for %d query phecode(s); excluded from controls: %s",
                    length(no_pos), paste(no_pos, collapse = ", ")))

  # Enumerate the valid candidate universe once (categories attached), then
  # explode it per query endpoint so each query has its own anchored pool.
  valid_universe <- BuildValidUniverse(query_phecodes, dt_nodes, dt_true)
  univ_anchored  <- AnchorByQueryEndpoint(valid_universe, query_phecodes)
  true_anchored  <- AnchorByQueryEndpoint(dt_true[, .(u, v)], query_phecodes)

  out <- vector("list", nrow(n_pos_tab))
  for (i in seq_len(nrow(n_pos_tab))) {
    q     <- n_pos_tab$phecode_query[i]
    n_pos <- n_pos_tab$n_pos[i]
    pool  <- univ_anchored[phecode_query == q]

    out[[i]] <- rbindlist(list(
      SampleQueryRatioCells(pool[, .(u, v)],
        n_pos, ratios, "random", q),
      SampleQueryRatioCells(pool[!is.na(cat_u) & !is.na(cat_v) & cat_u != cat_v, .(u, v)],
        n_pos, ratios, "non_adjacent", q),
      SampleQueryRatioCells(pool[!is.na(cat_u) & !is.na(cat_v) & cat_u == cat_v, .(u, v)],
        n_pos, ratios, "intracategory", q),
      SampleDegreeMatchedQuery(q, true_anchored[phecode_query == q],
        all_phecodes, dt_true, n_pos, ratios)
    ), use.names = TRUE, fill = TRUE)
  }
  return(rbindlist(out, use.names = TRUE, fill = TRUE))
}

# =========================================================================
# Operating-characteristics evaluation
# Append to functions.R. Requires: data.table, ggplot2.
# =========================================================================

# -------------------------------------------------------------------------
# Predicted scores: canonicalize BEFORE aggregating
# -------------------------------------------------------------------------
# embed_clinical pairs are (query, retrieved) i.e. directional and NOT
# canonical, while the validation edges are canonical (u < v). Canonicalizing
# after aggregation would require re-aggregating mean across the two
# directions; doing it first is exact and cheaper. If both directions exist
# for a pair, they collapse into one canonical key and the agg sees both rows.
AggregatePredicted <- function(embed_clinical,
                               score_col = "score_cosine",
                               u_col = "u", v_col = "v") {
  dt <- copy(embed_clinical)
  setnames(dt, c(u_col, v_col, score_col), c("u", "v", ".s"))
  # Preserve the query phecode (the u_col endpoint) BEFORE canonicalization
  # swaps u/v, so per-query identity survives into the aggregated scores. The
  # same canonical (u, v) retrieved from two different queries stays separate.
  dt[, phecode_query := u]
  dt[, `:=`(u = fifelse(u < v, u, v), v = fifelse(u < v, v, u))]
  dt <- dt[u != v]
  return(dt[, .(score_max  = max(.s,  na.rm = TRUE),
         score_mean = mean(.s, na.rm = TRUE),
         score_min  = min(.s,  na.rm = TRUE)),
     by = .(phecode_query, u, v)]
  )
}

# -------------------------------------------------------------------------
# Broadcast positives into every (stratum, grid_size) cell
# -------------------------------------------------------------------------
# Positives are the SAME true-edge set in every cell; negatives are
# cell-specific. This replicates the n_cases positives into each cell so each
# (stratum, grid_size) has a complete labelled set. Resolves the "true edges
# have NA grid_size" issue: the grid_size on a positive row is the cell it is
# being evaluated *in*, not a property of the edge.
BroadcastPositives <- function(dt_true, dt_synthetic_control) {
  pos <- copy(dt_true)
  pos[, `:=`(u = fifelse(u < v, u, v), v = fifelse(u < v, v, u))]
  pos <- unique(pos[u != v, .(u, v)])

  # Anchor positives per query endpoint (same convention as the negatives), so
  # each query's positive set is replicated into every (stratum, grid_size) cell
  # for that query. A query-query positive is counted under both queries.
  query_phecodes <- unique(dt_synthetic_control$phecode_query)
  pos_q <- AnchorByQueryEndpoint(pos, query_phecodes)
  pos_q[, true_edge := 1L]

  cells <- unique(dt_synthetic_control[, .(phecode_query, stratum, grid_size)])
  pos_b <- merge(cells, pos_q, by = "phecode_query", allow.cartesian = TRUE)[
    , .(phecode_query, u, v, stratum, grid_size, true_edge)]

  neg <- dt_synthetic_control[, .(phecode_query, u, v, stratum, grid_size, true_edge = 0L)]
  return(rbindlist(list(pos_b, neg), use.names = TRUE))
}

# -------------------------------------------------------------------------
# Attach scores, melt to long over score_type, handle unscored pairs
# -------------------------------------------------------------------------
# unscored = "floor": pair absent from dt_predicted is treated as a low score
#   (below every observed score). This is the correct retrieval-eval behaviour
#   when embed_clinical is top-k: "not retrieved" == "predicted negative".
#   The exact floor value is irrelevant to AUROC/AP as long as it sits below
#   all real scores; ties among floored pairs are handled by the rank/AP code.
# unscored = "drop": evaluate only scored pairs. Almost always WRONG here
#   (drops the bulk of negatives -> severe selection bias). Provided for the
#   dense-score case (embed_clinical is a full similarity matrix, no NAs).
ScoreEdges <- function(dt_eval, dt_predicted,
                       score_cols  = c("score_max", "score_mean", "score_min"),
                       unscored    = c("floor", "drop"),
                       floor_value = NULL) {
  unscored <- match.arg(unscored)

  # dt_eval already carries phecode_query: positives are anchored per query in
  # BroadcastPositives and negatives per query in GenerateSyntheticControls (a
  # pair with two query endpoints appears under both). Scores therefore attach
  # on the full per-query key; a pair a given query never retrieved stays
  # unscored and is floored below (from that query's view, a predicted negative).
  dt <- merge(dt_eval, dt_predicted,
              by = c("phecode_query", "u", "v"), all.x = TRUE)

  dt_long <- melt(
    dt,
    id.vars      = c("u", "v", "phecode_query", "stratum", "grid_size", "true_edge"),
    measure.vars = score_cols,
    variable.name = "score_type", value.name = "score"
  )
  dt_long[, score_type := sub("^score_", "", score_type)]

  if (unscored == "drop") {
    n_drop <- dt_long[is.na(score), .N]
    if (n_drop > 0L)
      warning(sprintf("ScoreEdges(drop): removing %d unscored edge-rows.", n_drop))
    dt_long <- dt_long[!is.na(score)]
  } else {
    dt_long[, score := {
      obs <- score[!is.na(score)]
      fv  <- if (!is.null(floor_value)) floor_value
             else if (length(obs)) min(obs) - 1 else 0
      fifelse(is.na(score), fv, score)
    }, by = score_type]
  }
  return(dt_long)
}

# -------------------------------------------------------------------------
# Scalar metrics
# -------------------------------------------------------------------------
# AUROC via the Mann-Whitney rank identity: exact, tie-corrected (midranks via
# frank average), and equals the trapezoidal ROC-AUC. Fast under data.table by.
Auroc <- function(score, label) {
  n_pos <- sum(label == 1L); n_neg <- sum(label == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  r <- frank(score, ties.method = "average")
  return((sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
}

# AUPRC as Average Precision (step-wise, no linear interpolation). This is the
# sklearn average_precision_score estimator and the one that is NOT optimistic
# the way linear PR interpolation is (Davis & Goadrich 2006). Ties handled by
# collapsing to distinct score thresholds. For a second opinion / canonical
# numbers, validate against precrec::evalmod or PRROC::pr.curve.
Auprc <- function(score, label) {
  n_pos <- sum(label == 1L)
  if (n_pos == 0L) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  s <- score[ord]; y <- label[ord]
  tp <- cumsum(y == 1L); fp <- cumsum(y == 0L)
  keep <- c(s[-length(s)] != s[-1L], TRUE)        # last index of each tie-run
  tp <- tp[keep]; fp <- fp[keep]
  precision <- tp / (tp + fp)
  recall    <- tp / n_pos
  d_recall  <- c(recall[1L], diff(recall))
  return(sum(d_recall * precision))
}

ComputeAUC <- function(dt_scored) {
  out <- dt_scored[, .(
    n_pos      = sum(true_edge == 1L),
    n_neg      = sum(true_edge == 0L),
    auroc      = Auroc(score, true_edge),
    auprc      = Auprc(score, true_edge)
  ), by = .(phecode_query, stratum, grid_size, score_type)]
  out[, `:=`(prevalence    = n_pos / (n_pos + n_neg),   # = AUPRC no-skill baseline
             ratio_realized = n_neg / n_pos)]           # may be < grid_size/n_cases
  return(out)                                           #   for saturated strata
}

# -------------------------------------------------------------------------
# Threshold sweep: sensitivity / specificity / PPV / NPV
# -------------------------------------------------------------------------
# Predict positive iff score >= threshold. Operating points at every distinct
# threshold within a group (the empirical ROC/PR support). Cumulative sort is
# O(m log m) per group. For very fine output, pass a coarser `thresholds`.
RocPoints <- function(score, label, thresholds = NULL) {
  P <- sum(label == 1L); N <- sum(label == 0L)
  ord <- order(score, decreasing = TRUE)
  s <- score[ord]; y <- label[ord]
  tp <- cumsum(y == 1L); fp <- cumsum(y == 0L)
  keep <- c(s[-length(s)] != s[-1L], TRUE)
  thr <- s[keep]; TP <- tp[keep]; FP <- fp[keep]
  FN <- P - TP; TN <- N - FP
  res <- data.table(
    threshold   = thr,
    tp = TP, fp = FP, fn = FN, tn = TN,
    sensitivity = TP / P,
    specificity = TN / N,
    ppv = fifelse(TP + FP > 0L, TP / (TP + FP), NA_real_),
    npv = fifelse(TN + FN > 0L, TN / (TN + FN), NA_real_)
  )
  if (!is.null(thresholds)) {
    # snap to nearest supplied threshold (predict-positive iff score >= t)
    idx <- findInterval(thresholds, rev(thr))   # rev -> ascending for findInterval
    res <- res[.N - idx + 1L]
    res[, threshold := thresholds]
  }
  return(res)
}

ComputeOperatingPoints <- function(dt_scored, thresholds = NULL) {
  return(dt_scored[, RocPoints(score, true_edge, thresholds),
            by = .(phecode_query, stratum, grid_size, score_type)])
}

# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------
# AUROC and AUPRC overlaid; facet stratum x score_type; grid_size (log) on x.
# Reference lines: 0.5 (AUROC no-skill) and the prevalence curve (AUPRC
# no-skill). The latter falls as grid_size grows, so a flat raw AUPRC vs
# grid_size actually means IMPROVING skill-over-baseline -- read AUPRC against
# its dashed reference, not in absolute terms.
PlotAUC <- function(dt_auc, x_log = TRUE) {
  # dt_auc is per query phecode; macro-average across phenotypes so each
  # stratum x score_type cell has a single AUROC and AUPRC line. NA-safe: a
  # query with no positives (e.g. COPD) yields NA per-query AUC and is skipped.
  agg <- dt_auc[, .(
    auroc      = mean(auroc, na.rm = TRUE),
    auprc      = mean(auprc, na.rm = TRUE),
    prevalence = mean(prevalence, na.rm = TRUE)
  ), by = .(stratum, grid_size, score_type)]
  d <- melt(agg,
            id.vars = c("stratum", "grid_size", "score_type", "prevalence"),
            measure.vars = c("auroc", "auprc"),
            variable.name = "metric", value.name = "value")
  p <- ggplot(d, aes(grid_size, value, colour = metric)) +
    geom_hline(yintercept = 0.5, linetype = 3, colour = "grey60") +
    geom_line(data = agg, aes(grid_size, prevalence),
              inherit.aes = FALSE, linetype = 2, colour = "grey50") +
    geom_line() + geom_point(size = 1) +
    facet_grid(stratum ~ score_type) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "negative : positive ratio (grid_size)",
         y = "metric (mean over phenotypes)", colour = NULL,
         caption = paste("lines = mean across query phenotypes;",
                         "dotted = AUROC no-skill (0.5);",
                         "dashed = AUPRC no-skill (mean prevalence)")) +
    theme_minimal(base_size = 11)
  if (x_log) p <- p + scale_x_log10()
  return(p)
}

# -------------------------------------------------------------------------
# Per-query colour palette & phenotype labels
# -------------------------------------------------------------------------
# The fixed per-query colour palette (QUERY_PHECODE_COLORS) and human-readable
# phenotype labels (QUERY_PHECODE_LABELS) used to live here as top-level vectors.
# They have been promoted to tracked pipeline targets (`query_phecode_colors`,
# `query_phecode_labels`) in `_targets.R` — see CLAUDE.md, which forbids
# top-level executable statements in R/ (tar_source() runs every file here).
# The plotting/table functions below now take `palette` / `labels` as explicit
# arguments, supplied from those targets.

# -------------------------------------------------------------------------
# Per-query AUC bar chart  (item 4)
# -------------------------------------------------------------------------
# One bar per query phecode, ordered by the chosen AUC metric (descending),
# coloured with the fixed QUERY_PHECODE_COLORS palette. Single panel.
# Defaults to the score_mean / random-stratum / 1:1 (smallest grid) condition.
# `metric` selects AUROC or AUPRC; report both per the project's metric rule.
PlotAUCByQuery <- function(dt_auc,
                           metric       = c("auroc", "auprc"),
                           score        = "mean",
                           stratum_keep = "random",
                           grid_keep    = NULL,
                           palette,
                           labels) {
  metric <- match.arg(metric)
  d <- dt_auc[stratum == stratum_keep & score_type == score & !is.na(phecode_query)]
  if (nrow(d) == 0L)
    stop("No rows after filtering; check stratum/score_type values in dt_auc.")
  # 1:1 case-to-negative ratio == the smallest grid (multiple 1 -> n_cases
  # negatives == n_cases positives, pooled).
  if (is.null(grid_keep)) grid_keep <- min(d$grid_size)
  d <- d[grid_size == grid_keep]
  d[, value := get(metric)]

  # x-axis and colours stay keyed on phecode internally (so the fixed palette
  # matches), but ticks are relabelled to phenotype names; legend dropped since
  # the axis already names each bar.
  p <- ggplot(d, aes(x = stats::reorder(phecode_query, -value),
                     y = value, fill = phecode_query)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = palette, drop = FALSE, guide = "none") +
    scale_x_discrete(labels = labels) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "query phenotype", y = toupper(metric),
         title = sprintf("%s by query phenotype", toupper(metric)),
         subtitle = sprintf("%s stratum, score_%s, %s:1 case:control",
                            stratum_keep, score, grid_keep)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (metric == "auroc")
    p <- p + geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50")
  return(p)
}

# -------------------------------------------------------------------------
# Per-query operating-characteristic curves  (item 5)
# -------------------------------------------------------------------------
# For each query phecode, a 2x2 panel of sensitivity / specificity / PPV / NPV
# vs. score threshold. Returns a NAMED LIST of ggplot objects (one per query
# phecode); iterate to print/save, e.g.
#   plots <- PlotOperatingPointsByQuery(dt_roc)
#   for (q in names(plots)) ggsave(sprintf("roc_%s.png", q), plots[[q]])
# Defaults to the score_mean / random-stratum / 1:1 (smallest grid) condition.
PlotOperatingPointsByQuery <- function(dt_roc,
                                       score        = "mean",
                                       stratum_keep = "random",
                                       grid_keep    = NULL,
                                       labels) {
  d <- dt_roc[stratum == stratum_keep & score_type == score & !is.na(phecode_query)]
  if (nrow(d) == 0L)
    stop("No rows after filtering; check stratum/score_type values in dt_roc.")
  if (is.null(grid_keep)) grid_keep <- min(d$grid_size)
  d <- d[grid_size == grid_keep]

  d_long <- melt(
    d,
    id.vars       = c("phecode_query", "threshold"),
    measure.vars  = c("sensitivity", "specificity", "ppv", "npv"),
    variable.name = "metric_type", value.name = "value"
  )
  metric_labels <- c(sensitivity = "Sensitivity", specificity = "Specificity",
                     ppv = "PPV", npv = "NPV")
  d_long[, metric_type := factor(metric_type,
                                 levels = names(metric_labels),
                                 labels = metric_labels)]

  queries <- sort(unique(d_long$phecode_query))
  plots <- lapply(queries, function(q) {
    nm    <- labels[q]
    title <- if (is.na(nm))
      sprintf("Operating characteristics — phecode %s", q)
    else
      sprintf("Operating characteristics — %s (phecode %s)", nm, q)
    ggplot(d_long[phecode_query == q], aes(x = threshold, y = value)) +
      geom_step(direction = "vh") +
      facet_wrap(~ metric_type, nrow = 2, ncol = 2) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "score threshold", y = "metric",
           title = title,
           subtitle = sprintf("%s stratum, score_%s, %s:1 case:control",
                              stratum_keep, score, grid_keep)) +
      theme_minimal(base_size = 11)
  })
  # Name the returned list by human-readable phenotype (e.g. "T2DM"), sanitized
  # for use as a filename, falling back to the raw phecode when no label exists.
  # Downstream targets write roc_query_<name>.png straight from these names, so
  # the abbreviation — not the phecode number — becomes the filename.
  keys <- vapply(queries, function(q) {
    nm   <- labels[q]
    base <- if (is.na(nm)) q else nm
    gsub("[^A-Za-z0-9._-]", "_", base)
  }, character(1))
  names(plots) <- make.unique(keys, sep = "_")
  return(plots)
}

# -------------------------------------------------------------------------
# Headline AUC table
# -------------------------------------------------------------------------
# The per-query AUROC/AUPRC summary at one condition (default: random stratum,
# score_mean, 1:1 ratio == smallest grid), ordered by AUROC descending and
# annotated with the readable phenotype name. Returned as a data.table for
# writing/printing.
HeadlineAucTable <- function(dt_auc,
                             score        = "mean",
                             stratum_keep = "random",
                             grid_keep    = NULL,
                             labels) {
  d <- dt_auc[stratum == stratum_keep & score_type == score & !is.na(phecode_query)]
  if (nrow(d) == 0L)
    stop("No rows after filtering; check stratum/score_type values in dt_auc.")
  if (is.null(grid_keep)) grid_keep <- min(d$grid_size)
  d <- d[grid_size == grid_keep]
  d[, phenotype := labels[phecode_query]]
  setorder(d, -auroc)
  return(d[, .(phenotype, phecode_query, n_pos, n_neg,
               auroc, auprc, prevalence, ratio_realized,
               stratum, grid_size, score_type)])
}

# -------------------------------------------------------------------------
# Output writer
# -------------------------------------------------------------------------
WriteEvalOutputs <- function(dt_auc, dt_roc, plot_auc, outdir = "output/dx") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  f_auc  <- file.path(outdir, "operating_characteristics_auc.csv")
  f_roc  <- file.path(outdir, "operating_characteristics_thresholds.csv")
  f_plot <- file.path(outdir, "auc_vs_gridsize.png")
  fwrite(dt_auc, f_auc)
  fwrite(dt_roc, f_roc)
  ggplot2::ggsave(f_plot, plot_auc, width = 11, height = 8, dpi = 300)
  c(auc = f_auc, roc = f_roc, plot = f_plot)
}