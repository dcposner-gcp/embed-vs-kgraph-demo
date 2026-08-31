# R/simulate.R -- synthetic inputs for the evaluation pipeline
#
# WHAT: side-effect-free functions that build every file `_targets.R` reads --
# the public Phecode 1.2 vocabulary (fetched from the PheWAS package mirror on
# GitHub, sha256-pinned), a synthetic ICD-9-CM / ICD-10-CM typing dictionary
# standing in for the OHDSI/Athena export, a KESER-shaped gold-standard graph
# (`dict.wiki.keep` in an .Rdata), and retriever-shaped ranked ICD retrievals
# (`<query>_Diagnosis.csv`) -- from ONE latent geometry with a known, tunable
# signal. The evaluator never learns the data is simulated: the files land in a
# directory laid out exactly like `data/`, and `tar_make()` runs unchanged.
# WHY: the real inputs (VA-trained KESER graph, the private embedding-model retrievals,
# OHDSI vocabulary) cannot be redistributed, so the public demo and the
# calibration checks (snr = 0 => AUROC ~ 0.5; oracle => 1) run on these.
# Generative model (see docs/simulator.md): nodes get
# comorbidity vectors V (category + root + mid + leaf offsets) that generate the
# gold KG; the retriever sees U = lambda*V + (1-lambda)*H + eps/snr where H is a
# hierarchy-only view (lambda = 0 is the hierarchy-only baseline). Driver:
# scripts/simulate_data.R. No top-level execution -- tar_source() loads this.

# ---- Public Phecode 1.2 source ---------------------------------------------

PhecodeSourceSpec <- function() {
  base <- paste0("https://raw.githubusercontent.com/PheWAS/PheWAS/",
                 "25bbb2f692bd7a7ba755921ec1b3753558d7d319/data/")
  list(
    map = list(
      file = "phecode_map.rda", object = "phecode_map",
      url = paste0(base, "phecode_map.rda"),
      sha256 = "53bb50e57a154a834a7862163ffacd754bc956f65df0efd5fb70d1d1d183ab15"),
    info = list(
      file = "pheinfo.rda", object = "pheinfo",
      url = paste0(base, "pheinfo.rda"),
      sha256 = "c4ba0f6ef21d07183a2ecbc155eec7278a7f2bcb2a7f86bbd1ce59274485fa71"),
    citation = paste(
      "Phecode 1.2 (ICD-9-CM map + ICD-10-CM beta map, flat/rolled) as shipped in",
      "the PheWAS R package, github.com/PheWAS/PheWAS (GPL-3), data/ at commit",
      "25bbb2f. Cite Denny JC et al. Nat Biotechnol 2013;31:1102-11 and",
      "Wu P et al. JMIR Med Inform 2019;7:e14325 (phewascatalog.org).")
  )
}

# Download (once) + verify + load the two tables. Returns list(map, info, spec):
#   map  : vocabulary_id (ICD9CM | ICD10CM), code, phecode   -- one row per ICD
#   info : phecode, description, groupnum, group             -- one row per phecode
FetchPhecodeSource <- function(cache_dir) {
  spec <- PhecodeSourceSpec()
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  load_one <- function(s) {
    path <- file.path(cache_dir, s$file)
    if (!file.exists(path)) {
      cli::cli_alert_info("Downloading {s$file} -> {cache_dir}")
      utils::download.file(s$url, path, mode = "wb", quiet = TRUE)
    }
    got <- digest::digest(file = path, algo = "sha256")
    if (!identical(got, s$sha256))
      stop(sprintf("sha256 mismatch for %s: got %s, expected %s", path, got, s$sha256))
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    as.data.table(e[[s$object]])
  }
  map  <- load_one(spec$map)[, .(vocabulary_id, code, phecode)]
  info <- load_one(spec$info)[, .(phecode, description, groupnum = as.integer(groupnum), group)]
  # The ICD-10-CM beta map carries ~120 chapter-block headings ("A00-A09") that
  # are not codes -- the pipeline's SplitRanges() would expand each into every
  # code of the block. Drop them. Also drop the handful of map rows pointing at
  # phecodes with no definition, so the node dictionary and the map agree.
  map <- map[!grepl("-", code, fixed = TRUE)]
  map <- unique(map[phecode %chin% info$phecode])
  list(map = map, info = unique(info), spec = spec)
}

# ---- Dotted hierarchy helpers (same convention as R/functions.R) -------------

# "250.22" -> "250.2" -> "250" -> NA
DottedParent <- function(id) {
  p <- ifelse(grepl(".", id, fixed = TRUE), substr(id, 1L, nchar(id) - 1L), NA_character_)
  sub("\\.$", "", p)
}

# Long table (node, ancestor) over every ancestor level.
DottedLineage <- function(ids) {
  ids <- unique(ids)
  out <- list()
  cur <- data.table(node = ids, ancestor = DottedParent(ids))
  cur <- cur[!is.na(ancestor)]
  while (nrow(cur) > 0L) {
    out[[length(out) + 1L]] <- copy(cur)
    cur[, ancestor := DottedParent(ancestor)]
    cur <- cur[!is.na(ancestor)]
  }
  if (length(out)) rbindlist(out) else data.table(node = character(), ancestor = character())
}

DottedDepth <- function(id) nchar(sub("^[^.]*\\.?", "", id))

# Mirrors IsAncestorDescendant() in R/functions.R (string-prefix rule).
IsDottedAncestor <- function(u, v) startsWith(v, u) & (nchar(v) > nchar(u))

IsHierarchyPair <- function(u, v) IsDottedAncestor(u, v) | IsDottedAncestor(v, u)

# Canonical undirected pairs: u < v (R string comparison, as the pipeline does),
# no self-loops, no duplicates.
CanonPairs <- function(dt) {
  dt <- copy(dt)
  dt[, `:=`(u = fifelse(u < v, u, v), v = fifelse(u < v, v, u))]
  unique(dt[u != v])
}

# ---- Vocabulary ---------------------------------------------------------------

BuildSimVocabulary <- function(src) {
  nodes <- src$info[, .(node_id = phecode, label = description,
                        category = group, category_number = groupnum)]
  setorder(nodes, node_id)
  nodes[, root := sub("\\..*", "", node_id)]
  nodes[, depth := DottedDepth(node_id)]
  par <- DottedParent(nodes$node_id)
  nodes[, parent_id := fifelse(par %chin% node_id, par, NA_character_)]
  n_icd <- src$map[, .(n_icd = .N), by = phecode]
  nodes[n_icd, n_icd := i.n_icd, on = .(node_id = phecode)]
  nodes[is.na(n_icd), n_icd := 0L]
  nodes[, is_parent := node_id %chin% DottedLineage(node_id)$ancestor]
  nodes[]
}

# ---- Latent geometry ----------------------------------------------------------

# One vector per node: category offset + one offset per lineage level (root,
# mid, leaf), each drawn ONCE per key so siblings share their ancestors' terms.
# Rows are L2-normalised, so tcrossprod() gives cosines.
SimulateLatent <- function(nodes, d = 32, sd_category = 1, sd_root = 1,
                           sd_mid = 0.7, sd_leaf = 0.5, seed = 1L) {
  set.seed(seed)
  n <- nrow(nodes)
  keyed <- function(keys, sd) {
    u <- unique(keys[!is.na(keys)])
    Z <- matrix(rnorm(length(u) * d, 0, sd), length(u), d, dimnames = list(u, NULL))
    out <- matrix(0, n, d)
    ok <- !is.na(keys)
    out[ok, ] <- Z[keys[ok], , drop = FALSE]
    out
  }
  M <- keyed(nodes$category, sd_category)
  lin <- rbind(DottedLineage(nodes$node_id),
               data.table(node = nodes$node_id, ancestor = nodes$node_id))
  lin[, adepth := DottedDepth(ancestor)]
  sds <- c(sd_root, sd_mid, sd_leaf)
  for (lv in 0:2) {
    key <- lin[adepth == lv][match(nodes$node_id, node), ancestor]
    M <- M + keyed(key, sds[lv + 1L])
  }
  M <- M / sqrt(rowSums(M^2))
  rownames(M) <- nodes$node_id
  M
}

# ---- Gold-standard graph -----------------------------------------------------

# Theoretical hierarchy closure of an edge set: for each edge (a, b), every
# (x, y) with x in lineage(a) and y in lineage(b); returns the ones NOT already
# present (canonical, hierarchy pairs excluded). Same construction the
# pipeline's CheckPhecodeHierarchy() uses to score closure compliance.
InheritedEdges <- function(edges, ids) {
  lin <- DottedLineage(ids)[ancestor %chin% ids]
  lin <- rbind(lin, data.table(node = ids, ancestor = ids))
  e <- edges[, .(id = .I, u, v)]
  a <- lin[e, on = .(node = u), .(id, x = ancestor), allow.cartesian = TRUE]
  b <- lin[e, on = .(node = v), .(id, y = ancestor), allow.cartesian = TRUE]
  th <- a[b, on = "id", allow.cartesian = TRUE][x != y, .(u = x, v = y)]
  th <- CanonPairs(th)
  th <- th[!IsHierarchyPair(u, v)]
  th[!edges, on = .(u, v)]
}

# Base edges from the latent geometry, then the flat-mapping quirk (only a
# fraction `closure` of hierarchy-inherited edges exist), parent dropout (some
# parent nodes are absent from the co-occurrence graph although their children
# are present -- what ExpandMissingParents() repairs), and uniform rewiring of
# `noise_frac` edges (edges no geometry explains).
#   edge_rule = "logistic": P(edge) = sigmoid(alpha + beta*cos + gamma*same_cat
#                          + pop_i + pop_j), alpha solved for mean_degree.
#   edge_rule = "threshold": edge iff cos + gamma*same_cat + pop_i + pop_j > tau
#                          (deterministic; with gamma = sd_pop = 0 the oracle
#                          retriever reaches AUROC = 1).
#   edge_rule = "loglinear": the Xu et al. generative model (SimulateLogLinearGraph):
#                          patient sequences -> co-occurrence -> PMI -> BH test. `pop`
#                          becomes log code prevalence; mean_degree/beta/gamma unused;
#                          the graph's density comes from n_patients, seq_len, kappa, fdr.
SimulateGoldEdges <- function(nodes, V, edge_rule = c("logistic", "threshold", "loglinear"),
                              mean_degree = 12, beta = 6, gamma = 1, sd_pop = 0.8,
                              closure = 0.34, parent_dropout = 0.3, noise_frac = 0.05,
                              protect = character(), seed = 2L, pop = NULL,
                              n_patients = 20000L, seq_len = 100L, kappa = 1.5, fdr = 0.05,
                              min_count = 5L) {
  edge_rule <- match.arg(edge_rule)
  set.seed(seed)
  ids <- nodes$node_id
  n <- length(ids)
  C <- tcrossprod(V)
  ij <- which(upper.tri(C), arr.ind = TRUE)
  i <- ij[, 1L]; j <- ij[, 2L]
  keep <- !IsHierarchyPair(ids[i], ids[j])
  i <- i[keep]; j <- j[keep]
  cosine <- C[cbind(i, j)]
  if (is.null(pop)) pop <- rnorm(n, 0, sd_pop)
  same <- as.numeric(nodes$category_number[i] == nodes$category_number[j])
  loglinear <- NULL
  if (edge_rule == "loglinear") {
    loglinear <- SimulateLogLinearGraph(V, logfreq = pop, n_patients = n_patients, seq_len = seq_len,
                                        kappa = kappa, fdr = fdr, min_count = min_count,
                                        seed = seed + 100L)
    hit <- NULL
    threshold <- NA_real_
  } else if (edge_rule == "logistic") {
    eta <- beta * cosine + gamma * same + pop[i] + pop[j]
    f <- function(a) mean(plogis(a + eta)) * (n - 1) - mean_degree
    alpha <- uniroot(f, c(-40, 20))$root
    hit <- runif(length(eta)) < plogis(alpha + eta)
    threshold <- alpha
  } else {
    eta <- cosine + gamma * same + pop[i] + pop[j]
    m <- round(mean_degree * n / 2)
    tau <- sort(eta, decreasing = TRUE)[m]
    hit <- eta >= tau
    threshold <- tau
  }
  base <- if (is.null(hit)) loglinear$edges[!IsHierarchyPair(u, v)] else
    CanonPairs(data.table(u = ids[i][hit], v = ids[j][hit]))
  base[, src := "base"]

  inherited <- InheritedEdges(base[, .(u, v)], ids)
  inherited <- inherited[runif(.N) < closure]
  inherited[, src := "closure"]
  edges <- rbind(base, inherited)

  parents <- nodes[is_parent == TRUE & !(node_id %chin% protect), node_id]
  dropped <- if (length(parents)) sample(parents, round(parent_dropout * length(parents))) else character()
  edges <- edges[!(u %chin% dropped | v %chin% dropped)]

  m_noise <- round(noise_frac * nrow(edges))
  if (m_noise > 0L) {
    edges <- edges[sample(.N, .N - m_noise)]
    rewired <- data.table(u = character(), v = character())
    while (nrow(rewired) < m_noise) {
      cand <- data.table(u = sample(ids, 2L * m_noise, replace = TRUE),
                         v = sample(ids, 2L * m_noise, replace = TRUE))
      cand <- CanonPairs(cand)[!IsHierarchyPair(u, v)]
      cand <- cand[!edges, on = .(u, v)][!rewired, on = .(u, v)]
      rewired <- unique(rbind(rewired, cand))
    }
    rewired <- rewired[seq_len(m_noise)]
    rewired[, src := "noise"]
    edges <- rbind(edges, rewired)
  }
  setorder(edges, u, v)
  list(edges = edges[], threshold = threshold, dropped_parents = dropped, pop = pop,
       loglinear = if (!is.null(loglinear)) loglinear[c("stats", "code_counts")] else NULL)
}

# Query-level guarantees: a designated zero-positive query (exercises the
# exclusion path), a minimum degree for every other query (top-up with the
# highest-cosine non-edges, i.e. consistent with the geometry), and forced
# query-query positives (exercise "counted under both queries").
ApplyQueryConstraints <- function(edges, V, query_nodes, zero_positive = character(),
                                  force_pairs = list(), min_degree = 5L) {
  edges <- copy(edges)
  edges <- edges[!(u %chin% zero_positive | v %chin% zero_positive)]
  ids <- rownames(V)
  for (q in setdiff(query_nodes, zero_positive)) {
    deg <- edges[u == q | v == q, .N]
    if (deg >= min_degree) next
    s <- as.vector(V[q, , drop = FALSE] %*% t(V))
    cand <- data.table(u = q, v = ids, s = s)[v != q]
    cand <- CanonPairs(cand)[!IsHierarchyPair(u, v)][!edges, on = .(u, v)]
    setorder(cand, -s)
    add <- cand[seq_len(min_degree - deg), .(u, v)]
    add[, src := "topup"]
    edges <- rbind(edges, add)
  }
  if (length(force_pairs)) {
    fp <- rbindlist(lapply(force_pairs, function(p) data.table(u = p[1L], v = p[2L])))
    fp <- CanonPairs(fp)[!(u %chin% zero_positive | v %chin% zero_positive)]
    fp <- fp[!edges, on = .(u, v)]
    if (nrow(fp)) { fp[, src := "forced"]; edges <- rbind(edges, fp) }
  }
  setorder(edges, u, v)
  edges[]
}

# ---- Retrieval -----------------------------------------------------------------

# Retriever's embedding. snr = 0 -> pure noise; snr = Inf -> no noise.
SimulateRetrieverEmbedding <- function(V, H, lambda = 0.7, snr = 2, seed = 3L) {
  set.seed(seed)
  n <- nrow(V); d <- ncol(V)
  E <- matrix(rnorm(n * d), n, d)
  E <- E / sqrt(rowSums(E^2))
  signal <- lambda * V + (1 - lambda) * H
  U <- if (snr == 0) E else if (is.infinite(snr)) signal else signal + E / snr
  U <- U / sqrt(rowSums(U^2))
  rownames(U) <- rownames(V)
  U
}

# Top-k phecode-level retrieval per query. `Q` is a matrix of query vectors with
# rownames = query_id; `self_node` (named by query_id, NA allowed) is excluded.
# k per query ~ round(k_mean * lognormal(0, k_sdlog)); k_mean = Inf -> dense.
# Score is a rank-preserving map of cosine into (0, 1], like a text-embedding
# cosine.
ScorePhecodeRetrieval <- function(U, Q, self_node, k_mean = 1200, k_sdlog = 0.3, seed = 4L) {
  set.seed(seed)
  S <- Q %*% t(U)
  ids <- rownames(U)
  n <- length(ids)
  out <- lapply(rownames(Q), function(q) {
    s <- (S[q, ] + 1) / 2
    keep <- rep(TRUE, n)
    if (!is.na(self_node[[q]])) keep[match(self_node[[q]], ids)] <- FALSE
    k <- if (is.infinite(k_mean)) n else round(k_mean * exp(rnorm(1L, 0, k_sdlog)))
    k <- min(k, sum(keep))
    o <- order(s, decreasing = TRUE)
    o <- o[keep[o]][seq_len(k)]
    data.table(query_id = q, candidate_id = ids[o], score = s[o], rank = seq_len(k))
  })
  rbindlist(out)
}

# ---- ICD layer -------------------------------------------------------------------

# Synthetic ICD typing dictionary (stands in for Athena's CONCEPT.csv). Every
# code in the Phecode map with its vocabulary, plus -- to reproduce the real
# cross-vocabulary ambiguity Athena exposes (ICD-9 V-codes reuse strings that
# are ICD-10-CM external-cause codes) -- an ICD10CM twin for a fraction of the
# ICD-9 V-roots. Returns list(dict, ambiguous_codes).
BuildSimIcdDictionary <- function(src_map, ambiguity_frac = 0.5, seed = 5L) {
  set.seed(seed)
  dict <- unique(src_map[, .(concept_code = code, vocabulary_id)])
  v9 <- dict[vocabulary_id == "ICD9CM" & startsWith(concept_code, "V"), concept_code]
  roots <- sort(unique(sub("\\..*", "", v9)))
  pick <- if (length(roots)) sample(roots, round(ambiguity_frac * length(roots))) else character()
  twin <- dict[vocabulary_id == "ICD9CM" & sub("\\..*", "", concept_code) %chin% pick,
               .(concept_code, vocabulary_id = "ICD10CM")]
  list(dict = unique(rbind(dict, twin)), ambiguous_codes = sort(unique(twin$concept_code)))
}

PerturbDescription <- function(x) {
  k <- sample(3L, length(x), replace = TRUE)
  fifelse(k == 1L, paste0(x, ", unspecified"),
  fifelse(k == 2L, sub("\\s+\\S+$", "", x), toupper(x)))
}

# Map phecode-level retrievals back to retriever-shaped raw ICD rows with the
# quirks the harmonisation stage handles: several ICD codes per phecode, comma
# lists, hyphen ranges, cross-vocabulary ambiguous codes, unmapped codes, a
# rank-1 self code with score 1, and one "ICD-10-CM" garbage row. Phecodes with
# no eligible ICD code (possible under icd_unique_only) are silently dropped.
# `queries`: query_id, label, node_id, self_code. Returns one row per raw line
# with `kind` (code|list|range|ambiguous|unmapped|self|garbage) for the manifest.
SimulateIcdLayer <- function(pred, queries, src, ambiguous_codes,
                             lambda_extra = 0.8, sd_icd = 0.01, p_list = 0.04,
                             p_range = 0.01, p_ambiguous = 0.02, p_unmapped = 0.05,
                             p_desc_perturb = 0.3, icd_unique_only = FALSE, seed = 6L) {
  set.seed(seed)
  desc_of <- setNames(src$info$description, src$info$phecode)
  # The ICD-10-CM map is many-to-many: ~4% of codes map to several phecodes, so a
  # retrieved code can carry its score onto a phecode the retriever never ranked
  # (the harmonisation stage keeps one of them). icd_unique_only = TRUE samples
  # only codes with a single phecode and skips cross-vocabulary ambiguous codes,
  # which is what lets the oracle preset reach AUROC = 1 exactly.
  pool <- src$map
  if (icd_unique_only) {
    n_phe <- pool[, .(n_phe = uniqueN(phecode)), by = code]
    pool <- pool[code %chin% n_phe[n_phe == 1L, code] & !(code %chin% ambiguous_codes)]
  }
  codes_by_phe <- split(pool$code, pool$phecode)
  vocab_of <- setNames(src$map$vocabulary_id, src$map$code)
  phe_of_code <- setNames(src$map$phecode, src$map$code)

  # 1. expand each retrieved phecode into 1 + Pois(lambda_extra) of its ICD codes
  n_avail <- lengths(codes_by_phe[pred$candidate_id])
  n_take <- pmin(1L + rpois(nrow(pred), lambda_extra), n_avail)
  picked <- mapply(function(ph, k) {
    cc <- codes_by_phe[[ph]]
    if (length(cc) <= k) cc else sample(cc, k)
  }, pred$candidate_id, n_take, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  rows <- data.table(query_id = rep(pred$query_id, n_take),
                     phecode  = rep(pred$candidate_id, n_take),
                     code     = unlist(picked),
                     score    = rep(pred$score, n_take))
  rows[, score := pmin(pmax(score + rnorm(.N, 0, sd_icd), 1e-4), 0.9999)]
  rows[, description := tolower(desc_of[phecode])]
  rows[runif(.N) < p_desc_perturb, description := PerturbDescription(description)]
  rows[, kind := "code"]

  # 2. comma lists: collapse some multi-code phecodes into one "A,B" row
  grp <- rows[, .N, by = .(query_id, phecode)][N >= 2L]
  grp <- grp[runif(.N) < p_list, .(query_id, phecode)]
  if (nrow(grp)) {
    merged <- rows[grp, on = .(query_id, phecode)][
      , .(code = paste(code, collapse = ","), score = max(score),
          description = description[1L], kind = "list"),
      by = .(query_id, phecode)]
    rows <- rbind(rows[!grp, on = .(query_id, phecode)], merged, use.names = TRUE)
  }

  # 3. hyphen ranges over adjacent roots of the same vocabulary ("410-414", "I20-I25")
  roots_by_vocab <- lapply(split(src$map$code, src$map$vocabulary_id),
                           function(cc) sort(unique(sub("\\..*", "", cc))))
  cand <- which(rows$kind == "code" & runif(nrow(rows)) < p_range)
  if (length(cand)) {
    r1 <- sub("\\..*", "", rows$code[cand])
    vv <- vocab_of[rows$code[cand]]
    r2 <- mapply(function(r, v) {
      rs <- roots_by_vocab[[v]]
      rs <- rs[nchar(rs) == nchar(r)]
      i <- match(r, rs); span <- sample(3L, 1L)
      rs[if (i + span <= length(rs)) i + span else max(1L, i - span)]
    }, r1, vv, USE.NAMES = FALSE)
    lo <- pmin(r1, r2); hi <- pmax(r1, r2)
    ok <- lo != hi
    rows[cand[ok], `:=`(code = paste0(lo[ok], "-", hi[ok]),
                        description = paste("codes", lo[ok], "to", hi[ok]),
                        kind = "range")]
  }

  # 4. ambiguous (both vocabularies) and unmapped extras, proportional to list length
  synthetic_unmapped <- unique(c(
    "U07.1", "U09.9",
    paste0(sample(LETTERS, 80, TRUE), sprintf("%02d", sample(0:99, 80, TRUE)), ".",
           sample(0:9, 80, TRUE)),
    sprintf("%02d.%02d", sample(0:99, 20, TRUE), sample(0:99, 20, TRUE))))
  unmapped_pool <- setdiff(synthetic_unmapped, src$map$code)
  per_q <- rows[, .N, by = query_id]
  extras <- per_q[, {
    n_amb <- rbinom(1L, N, p_ambiguous)
    n_unm <- rbinom(1L, N, p_unmapped)
    amb <- if (n_amb > 0L && length(ambiguous_codes)) sample(ambiguous_codes, n_amb, TRUE) else character()
    unm <- if (n_unm > 0L) sample(unmapped_pool, n_unm, TRUE) else character()
    .(code = c(amb, unm),
      kind = c(rep("ambiguous", length(amb)), rep("unmapped", length(unm))))
  }, by = query_id]
  extras <- unique(extras)
  if (nrow(extras)) {
    extras[, phecode := fifelse(kind == "ambiguous", phe_of_code[code], NA_character_)]
    extras[kind == "ambiguous", `:=`(description = tolower(desc_of[phecode]),
                                     score = runif(.N, 0.35, 0.75))]
    extras[kind == "ambiguous" & runif(.N) < 0.5, description := PerturbDescription(description)]
    extras[kind == "unmapped", `:=`(description = "unmapped code", score = runif(.N, 0.2, 0.6))]
    rows <- rbind(rows, extras, use.names = TRUE, fill = TRUE)
  }

  # 5. fixed rows per query: rank-1 self code (score 1), one "ICD-10-CM" garbage
  #    row, and one unmapped range ("U07-U09", the provisional COVID block: it
  #    expands to nothing). The last one guarantees at least one hyphenated code
  #    per data set -- the pipeline's SplitRanges() assumes ranges exist.
  self <- queries[, .(query_id, phecode = node_id, code = self_code,
                      score = 1, description = tolower(label), kind = "self")]
  garbage <- queries[, .(query_id, phecode = NA_character_, code = "ICD-10-CM",
                         score = 0.5, description = "ICD-10-CM", kind = "garbage")]
  urange <- queries[, .(query_id, phecode = NA_character_, code = "U07-U09",
                        score = 0.3, description = "provisional codes u07-u09", kind = "unmapped")]
  rows <- rbind(rows, self, garbage, urange, use.names = TRUE, fill = TRUE)
  setorder(rows, query_id, -score)
  rows[, rank := seq_len(.N), by = query_id]
  rows[]
}

# ---- Writers (the exact files `_targets.R` reads) ------------------------------------

# dim/Phecode_map_v1_2_icd9_icd10cm.csv + dim/phecode_definitions1.2.csv.
# ICDString is the PHECODE description (the public map ships no ICD strings;
# ICD-10-CM descriptions are not redistributed) -- the pipeline only uses it in
# the ambiguity audit, where the simulated retrieval descriptions match it.
WriteSimPhecodeFiles <- function(src, out_dir) {
  dim_dir <- file.path(out_dir, "dim")
  dir.create(dim_dir, showWarnings = FALSE, recursive = TRUE)
  phemap <- src$map[src$info, on = "phecode", nomatch = NULL,
                    .(ICD = code, Phecode = phecode, ICDString = description,
                      PhecodeString = description, PhecodeCategory = group)]
  setorder(phemap, ICD, Phecode)
  f_map <- file.path(dim_dir, "Phecode_map_v1_2_icd9_icd10cm.csv")
  fwrite(phemap, f_map)
  defs <- src$info[, .(phecode, phenotype = description, category_number = groupnum, category = group)]
  setorder(defs, phecode)
  f_def <- file.path(dim_dir, "phecode_definitions1.2.csv")
  fwrite(defs, f_def)
  c(f_map, f_def)
}

# dim/icd_athena_2025/CONCEPT.csv -- Athena column layout, tab-separated,
# unquoted, but SYNTHETIC: derived from the public Phecode map, not OHDSI.
WriteSimConceptCsv <- function(dict, out_dir) {
  d <- file.path(out_dir, "dim", "icd_athena_2025")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  concept <- dict[, .(
    concept_id = seq_len(.N), concept_name = "synthetic",
    domain_id = "Condition", vocabulary_id,
    concept_class_id = fifelse(grepl(".", concept_code, fixed = TRUE), "n-char billing code", "3-char nonbill code"),
    standard_concept = "", concept_code,
    valid_start_date = "19700101", valid_end_date = "20991231", invalid_reason = "")]
  f <- file.path(d, "CONCEPT.csv")
  fwrite(concept, f, sep = "\t", quote = FALSE)
  writeLines(c(
    "SYNTHETIC -- not an OHDSI/Athena export.",
    "Built by scripts/simulate_data.R from the public Phecode 1.2 map so the",
    "pipeline's ICD typing stage (PopulateDictionary / GetCodeType) runs without",
    "a licensed vocabulary download. Only concept_code + vocabulary_id are real;",
    "every other column is a placeholder."),
    file.path(d, "README.txt"))
  f
}

# validation/dict-disease-disease-wiki-bin-eval.Rdata -- object `dict.wiki.keep`
# with PheCode:-prefixed phecode1 / phecode2, both directions (a neighbourhood
# dictionary), as LoadKeserDxDx() expects.
WriteSimGoldRdata <- function(edges, out_dir) {
  d <- file.path(out_dir, "validation")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  both <- rbind(edges[, .(a = u, b = v)], edges[, .(a = v, b = u)])
  setorder(both, a, b)
  dict.wiki.keep <- data.frame(phecode1 = paste0("PheCode:", both$a),
                               phecode2 = paste0("PheCode:", both$b),
                               stringsAsFactors = FALSE)
  f <- file.path(d, "dict-disease-disease-wiki-bin-eval.Rdata")
  save(dict.wiki.keep, file = f, version = 2)
  f
}

# prediction/<query_id>_Diagnosis.csv: rank, code, description, score_cosine.
WriteSimPredictions <- function(icd_rows, out_dir) {
  d <- file.path(out_dir, "prediction")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  vapply(sort(unique(icd_rows$query_id)), function(q) {
    f <- file.path(d, paste0(q, "_Diagnosis.csv"))
    fwrite(icd_rows[query_id == q, .(rank, code, description, score_cosine = round(score, 7))], f)
    f
  }, character(1L))
}

# SPEC §1 tables (docs/SPEC-inputs.md): nodes,
# gold_edges, queries, predictions (phecode level, i.e. post-harmonisation) and,
# when `U` is given, embeddings (the retriever's vectors, one wide row per node:
# node_id, dim_1 .. dim_d -- SPEC 1.5 flattened for CSV).
WriteSimSpecTables <- function(nodes, edges, queries, pred, out_dir, U = NULL) {
  d <- file.path(out_dir, "spec")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(U)) {
    emb <- as.data.table(round(U, 6))
    setnames(emb, paste0("dim_", seq_len(ncol(U))))
    emb[, node_id := rownames(U)]
    setcolorder(emb, "node_id")
    fwrite(emb, file.path(d, "embeddings.csv"))
  }
  fwrite(nodes[, .(node_id, label, category, parent_id, node_type = "dx")],
         file.path(d, "nodes.csv"))
  fwrite(edges[, .(u, v, source = src)], file.path(d, "gold_edges.csv"))
  fwrite(queries[, .(query_id, label, node_id)], file.path(d, "queries.csv"))
  fwrite(pred[, .(query_id, candidate_id, rank, score = round(score, 7))],
         file.path(d, "predictions.csv"))
  d
}

# ---- Log-linear (Xu et al.) gold graph ----------------------------------------------
# The dynamic log-linear topic model of Xu, Gan, Zhou, Shen, Lu & Cai (JRSS-B 88(2), 2026;
# github.com/junwei-lu/WordVec_Inference, simulation/fig1and5/SIM_THMS.R), which is
# the theory under KESER/ARCH-style co-occurrence graphs. Direction is latent
# vectors -> patient event sequences -> windowed co-occurrence -> SPPMI -> test:
#   discourse  c_{t+1} = sqrt(a) c_t + sqrt(1-a) r_t,  r_t ~ N(0, I_p/p),  a = 1 - log(d)/p^2
#   emission   w_t ~ Multinomial(softmax(V c_t))
#   counts     C_ij = # (w_t, w_{t+u}) pairs with 1 <= u <= q (symmetric; diagonal x2)
#   PMI_ij     = log(max(eta, c * C_ij / (r_i r_j))),  c = (2qT - q^2 - q) n,  r = rowSums(C)
#   identity   PMI ~ alpha_p (V - mu)(V - mu)^T   (their eq. 2.5), so cosine of the
#              rank-p SVD embedding estimates cosine of the centred latent vectors
#   edge       PMI_ij / sd_ij significantly > 0 after Benjamini-Hochberg over all pairs
#              (their construct_KG.R), with sd_ij the patient-level empirical
#              estimator of Theorem 3.6 (SIM_THMS.R, "PMI_entry_sd_co"):
#              X_pat = C_pat,ij/mean(C_ij) - r_pat,i/mean(r_i) - r_pat,j/mean(r_j),
#              sd_ij = sd_pat(X) / sqrt(n).
# Differences from their script, on purpose: the walk starts at its stationary
# distribution (their code starts at 0 and burns in); V's direction comes from the
# hierarchy-structured latent vectors of this simulator and its row norms carry
# code prevalence -- under the model the marginal frequency of code w is
# proportional to exp(|V_w|^2 / (2p)), so |V_w|^2 = p*kappa^2 + 2p*logfreq_w gives
# frequency proportional to exp(logfreq_w) and a logit scale of about kappa.
# In this model the graph is sparse because of finite-sample POWER (every pair has a
# nonzero true PMI); rare codes get fewer edges. That is the property the logistic
# rule cannot reproduce.
SimulateLogLinearGraph <- function(V, logfreq, n_patients = 20000L, seq_len = 100L, q = 5L,
                                   kappa = 1, fdr = 0.05, min_count = 5L, eta = 1e-3, seed = 7L,
                                   chunk = 100L, verbose = TRUE) {
  set.seed(seed)
  d <- nrow(V); p <- ncol(V); Tn <- as.integer(seq_len); q <- as.integer(q)
  ids <- rownames(V)
  a <- 1 - log(d) / p^2
  norm2 <- pmax(p * kappa^2 + 2 * p * logfreq, 0.01)
  Vll <- V * sqrt(norm2)                       # rows of V are unit length
  VllT <- t(Vll)
  n_patients <- as.integer(n_patients)
  n_chunks <- ceiling(n_patients / chunk)
  pair_parts <- vector("list", n_chunks)
  code_counts <- integer(d)
  t0 <- Sys.time()
  for (ch in seq_len(n_chunks)) {
    n_c <- min(chunk, n_patients - (ch - 1L) * chunk)
    # stationary AR(1) discourse walk, all n_c patients at once: T steps of (n_c x p)
    Z <- array(0, c(n_c, p, Tn))
    z <- matrix(rnorm(n_c * p, 0, 1 / sqrt(p)), n_c, p)
    for (t in seq_len(Tn)) {
      if (t > 1L) z <- sqrt(a) * z + sqrt(1 - a) * matrix(rnorm(n_c * p, 0, 1 / sqrt(p)), n_c, p)
      Z[, , t] <- z
    }
    Zm <- matrix(aperm(Z, c(1, 3, 2)), n_c * Tn, p)          # row = (patient, t)
    L <- Zm %*% VllT                                           # logits, (n_c*T) x d
    L <- exp(L - matrixStats_rowMax(L))
    P <- L / rowSums(L)
    cum <- t(apply(P, 1L, cumsum))
    w <- rowSums(cum < runif(nrow(cum))) + 1L
    w[w > d] <- d
    W <- matrix(w, n_c, Tn)                                    # patient x time
    code_counts <- code_counts + tabulate(w, d)
    pat <- (ch - 1L) * chunk + seq_len(n_c)
    prs <- rbindlist(lapply(seq_len(q), function(u) {
      A <- W[, seq_len(Tn - u), drop = FALSE]; B <- W[, u + seq_len(Tn - u), drop = FALSE]
      data.table(pat = rep(pat, Tn - u), a = pmin(as.vector(A), as.vector(B)),
                 b = pmax(as.vector(A), as.vector(B)))
    }))
    pair_parts[[ch]] <- prs[, .(n = .N), by = .(pat, a, b)]
    if (verbose && (ch %% 20L == 0L || ch == n_chunks))
      cli::cli_alert_info("log-linear cohort: {ch}/{n_chunks} chunks, {round(as.numeric(Sys.time() - t0, units = 'secs'))} s")
  }
  PP <- rbindlist(pair_parts); rm(pair_parts)
  PP[a == b, n := 2L * n]                                      # their coocur_cal counts self-pairs twice
  Cdt <- PP[, .(C = sum(n)), by = .(a, b)]
  # dense symmetric co-occurrence matrix, row sums, total
  C <- matrix(0, d, d)
  C[cbind(Cdt$a, Cdt$b)] <- Cdt$C
  C[cbind(Cdt$b, Cdt$a)] <- Cdt$C
  r <- rowSums(C)
  ctot <- (2 * q * Tn - q^2 - q) * n_patients
  # empirical PMI on the tested pairs (i < j, C_ij > 0); everything else is floored
  tested <- Cdt[a < b]
  tested[, pmi := log(pmax(eta, ctot * C / (r[a] * r[b])))]
  # Theorem 3.6 patient-level variance of the empirical PMI entry
  rp <- rbind(PP[a != b, .(pat, i = a, n)], PP[a != b, .(pat, i = b, n)], PP[a == b, .(pat, i = a, n)])
  rp <- rp[, .(rn = sum(n)), by = .(pat, i)]
  rp[, B := rn / (r[i] / n_patients)]                          # r_pat,i / mean(r_i)
  Bm <- Matrix::sparseMatrix(i = rp$pat, j = rp$i, x = rp$B, dims = c(n_patients, d))
  G <- as.matrix(Matrix::crossprod(Bm)) / n_patients            # E[B_i B_j]
  PP[tested, Cij := i.C, on = .(a, b)]
  PP <- PP[a < b]
  PP[, A := n / (Cij / n_patients)]                             # C_pat,ij / mean(C_ij)
  PP[rp, Bi := i.B, on = .(pat, a = i)]
  PP[rp, Bj := i.B, on = .(pat, b = i)]
  mom <- PP[, .(EA2 = sum(A^2) / n_patients, EABi = sum(A * Bi) / n_patients,
                EABj = sum(A * Bj) / n_patients), by = .(a, b)]
  tested[mom, `:=`(EA2 = i.EA2, EABi = i.EABi, EABj = i.EABj), on = .(a, b)]
  tested[, EX2 := EA2 + G[cbind(a, a)] + G[cbind(b, b)] - 2 * EABi - 2 * EABj + 2 * G[cbind(a, b)]]
  tested[, sd := sqrt(pmax(EX2 - 1, 0) * n_patients / (n_patients - 1) / n_patients)]
  tested[, z := fifelse(sd > 0, pmi / sd, NA_real_)]
  tested[, p_one := pnorm(-z)]
  M <- d * (d - 1) / 2
  # low-rank estimator and its cosine (their construct_KG.R), for diagnostics
  alpha_p <- 2 * sum((Tn - seq_len(q)) * a^seq_len(q)) / (p * (2 * q * Tn - q^2 - q))
  PMI <- matrix(log(eta), d, d)
  PMI[cbind(tested$a, tested$b)] <- tested$pmi; PMI[cbind(tested$b, tested$a)] <- tested$pmi
  diag(PMI) <- log(pmax(eta, fifelse(r > 0, ctot * diag(C) / r^2, eta)))   # codes never seen: floor
  eg <- eigen(PMI, symmetric = TRUE)
  Vhat <- sqrt(alpha_p) * eg$vectors[, seq_len(p)] %*% diag(sqrt(pmax(eg$values[seq_len(p)], 0)))
  Vhat_n <- Vhat / pmax(sqrt(rowSums(Vhat^2)), 1e-12)
  tested[, cos_hat := rowSums(Vhat_n[a, ] * Vhat_n[b, ])]
  tested[, `:=`(u = ids[a], v = ids[b])]
  out <- tested[, .(u, v, C, pmi, sd, z, p_one, cos_hat)]
  edges <- LogLinearEdges(out, min_count = min_count, fdr = fdr, M = M)
  out[edges, edge := TRUE, on = .(u, v)]; out[is.na(edge), edge := FALSE]
  list(edges = edges, tested = out, code_counts = setNames(code_counts, ids),
       Vhat = `rownames<-`(Vhat, ids), alpha_p = alpha_p,
       stats = list(n_patients = n_patients, seq_len = Tn, q = q, kappa = kappa, fdr = fdr,
                    min_count = min_count, n_events = n_patients * Tn, n_tested = nrow(tested),
                    n_eligible = out[C >= min_count, .N], n_edges = nrow(edges),
                    frac_pairs_observed = nrow(tested) / M, alpha = a, alpha_p = alpha_p,
                    seconds = round(as.numeric(Sys.time() - t0, units = "secs"))))
}

# Edge call from the tested-pairs table: one-sided z-test on the empirical PMI,
# Benjamini-Hochberg at `fdr` over ALL d(d-1)/2 pairs (their construct_KG.R;
# unobserved pairs count as p = 1). `min_count` is ours: the normal approximation
# of Theorem 3.6 needs a non-trivial count, and without it a single co-occurrence
# of two rare codes is called significant (z ~ 5 at C = 1).
LogLinearEdges <- function(tested, min_count = 5L, fdr = 0.05, M) {
  t <- tested[C >= min_count & !is.na(z), .(u, v, pmi, p_one)]
  setorder(t, p_one)
  t[, k := .I]
  idx <- t[p_one <= k / M * fdr, if (.N) max(k) else 0L]
  CanonPairs(t[k <= idx & pmi > 0, .(u, v)])
}

# row-wise max in base R (max.col is C code; ties are irrelevant here)
matrixStats_rowMax <- function(m) m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
