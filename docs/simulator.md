# Simulator

Generates a complete `data/`-shaped input tree — plus the phecode-level spec
tables (`docs/SPEC-inputs.md`) — from Xu et al's statistical knowledge-graph model (PMID: 41982780) tuned to emulate Hong et al's KESER network (PMID: 34707226).
Files written to `data_sim/<preset>/` feed directly to `tar_make()`.

R functions: `R/simulate.R`
Driver script: `scripts/simulate_data.R` 
Driver options: presets, CLI, per-step seeds `seed+1 … seed+7`, manifest 
Run: `Rscript scripts/simulate_data.R --preset phecode`
(`--set key=value` overrides any preset parameter; default seed 20260829).

## 1. Design: two views of one vocabulary

The gold standard and the embedding are different *views* of the same latent
structure, not two samples of one process:

- `V` — latent **comorbidity** vectors, one per node. They generate the
  co-occurrence process and hence the gold KG.
- `H` — **hierarchy / lexical** vectors, the same construction with different
  variance shares (small leaf noise, larger category share, so siblings are
  near-identical). They stand in for what a description embedder sees.
- The retriever's embedding is `U = λ·V + (1 − λ)·H + ε / snr`
  (`SimulateRetrieverEmbedding()`). `λ` is how much co-occurrence structure
  the model learned beyond lexical similarity; `snr` is estimation noise.
  `λ = 0` is the **hierarchy-only baseline**; `snr = ∞, λ = 1` is the oracle.

## 2. Generative steps

1. **Vocabulary** — `FetchPhecodeSource()` + `BuildSimVocabulary()`. The real
   public Phecode 1.2 dictionary (~1.9k nodes): map + definitions downloaded
   from the PheWAS package mirror (cached in `data_sim/_source/`). 
   Dotted ids (`411`, `411.2`) give the hierarchy;
   categories come from the phecode groups. There is no synthetic-vocabulary
   mode — every preset runs on the real dictionary.
2. **Latent vectors** — `SimulateLatent(nodes, d, sd_category, sd_root,
   sd_mid, sd_leaf, seed)`: `V_i = μ_category + δ_root + δ_mid + ε_leaf`,
   each offset drawn once per key so siblings share their ancestors' terms;
   rows L2-normalised so `tcrossprod()` gives cosines. `H` uses the same
   function with the `h_sd_*` parameters.
3. **Gold KG** — 
  `SimulateGoldEdges()`, three `edge_rule` options with one output schema:
   - `logistic` (default): `P(edge_ij) = σ(α + β·cos(V_i,V_j) +
     γ·[same category] + pop_i + pop_j)`, `α` set to expected mean degree 
     hits `mean_degree`; `pop ~ N(0, sd_pop)` gives heavier-tailed degrees.
   - `threshold`: edge iff `cos + γ·[same cat] + pop_i + pop_j > τ`, `τ` set
     by `mean_degree`. Deterministic given `V` — with `γ = sd_pop = 0` 
     makes the oracle test exact; under the logistic rule the graph
     is stochastic given `V`, so even a perfect retriever AUROC < 1.
   - `loglinear`: the dynamic log-linear topic model of Xu, Gan, Zhou, Shen,
     Lu & Cai (JRSS-B 88(2), 2026; arXiv:2312.15611; code
     github.com/junwei-lu/WordVec_Inference). Implemented in full in
     `SimulateLogLinearGraph()`: per patient an AR(1) discourse
     `c_{t+1} = √a·c_t + √(1−a)·r_t`, events `w_t ~ Multinomial(softmax(V c_t))`,
     windowed co-occurrence (q = 5), shifted PMI, and an entrywise one-sided
     z-test (patient-level Theorem-3.6 variance) with BH over all pairs, plus
     a `min_count` eligibility filter. `mean_degree`/`beta`/`gamma` are
     unused; density is emergent from `n_patients`, `seq_len`, `kappa`,
     `fdr`. Sparsity here comes from finite-sample power — rare codes get
     fewer edges — which is the property the logistic rule cannot reproduce.

   Then simulate the **flat-mapping quirks** observed in the KESER knowledge-graph: 
   a) keep each hierarchy-inherited edge with probability `closure` (= 0.1 is ~40 % compliance under the pipeline's `CheckPhecodeHierarchy()`; the real KESER graph shows ~34 %); 
   b) drop a `parent_dropout` fraction of parent nodes from the graph that the
   pipeline's `ExpandMissingParents()` repairs; 
   query nodes are protected;
   c) rewire `noise_frac` of edges uniformly (edges no geometry explains, so the
   oracle needs `noise_frac = 0`). Every edge carries a `source` tag:
   `base`, `closure`, `noise`, `topup`, `forced`.
4. **Queries** — fixed, not sampled: the 11-query set the pipeline hard-codes
   (`QueryTable()` in the driver; CAD, COPD, COVID, Depression, Hearing loss,
   HF, Long COVID, MI, MS, RA, T2DM). COVID and Long COVID have no Phecode
   1.2 code: their retrieval is anchored on a stand-in node and their rank-1
   self code is an unmapped ICD-10 code (exercises the unmapped path).
   `ApplyQueryConstraints()` then guarantees: one designated **zero-positive**
   query (COPD — exercises the exclusion path), a minimum gold degree for
   every other query (topped up with the highest-cosine non-edges, i.e.
   consistent with the geometry), and forced **query–query gold pairs**
   (MI–CAD, T2DM–CAD — exercise "counted under both queries").
5. **Retrieval** — `ScorePhecodeRetrieval(U, Q, self_node, k_mean, k_sdlog)`:
   scores `cos(U_q, U_j)` mapped rank-preservingly into (0, 1], self code
   excluded; keep top-`k` per query with
   `k ~ round(k_mean · lognormal(0, k_sdlog))`, so recall@k varies and the
   floor path is exercised; `k_mean = Inf` ⇒ dense.
6. **ICD layer** — `SimulateIcdLayer()` maps retrieved phecodes back to fine
   ICD codes with the quirks the pipeline's harmonisation stage must handle:
   comma lists, dash ranges, cross-vocabulary ambiguous codes, unmapped
   codes, perturbed descriptions. Always on; the `phecode-oracle` preset
   neutralises the quirks. Note the oracle is exact only at the **phecode
   level** — through the raw-ICD layer it is not, for reasons that are
   properties of the harmonisation, not the evaluator: the Phecode 1.2
   ICD-10-CM map is many-to-many (~4 % of codes), cross-vocabulary ambiguous
   codes are dropped, and `ExpandMissingParents()` *adds* inherited positives
   for zero-degree parents. `phecode-oracle` neutralises the first two
   (`icd_unique_only`, no ambiguity) and reaches 1.000 in every stratum
   except `intracategory` (0.998) — the remaining gap is the third item.

## 3. Logistic vs log-linear gold graphs

Measured differences on identical `V` and prevalence at matched density
(20k patients × 100 events, κ = 1.5): degree tracks prevalence (Spearman
0.97 vs 0.81), 47 % vs 76 % same-category edges, max degree 152 vs 253,
transitivity 0.21 vs 0.26, edge Jaccard 0.20; with the same retriever the
log-linear graph is harder (mean AUROC 0.61 vs 0.68 at that density; 0.68 vs
0.73 at the KESER-like density of the `phecode-xu` preset, 5k patients,
mean degree ~14).

## 4. Presets

Defined in `Presets()` in `scripts/simulate_data.R` (not a config file). All
run on the real Phecode 1.2 vocabulary with the 11-query set; shared defaults
`d = 32, λ = 0.7, snr = 1.5, k_mean = 1200, mean_degree = 12`.

| preset | gold graph | deviations from `phecode` | expectation / use |
|---|---|---|---|
| `phecode` | logistic | — | README demo; realistic mid-range AUROC/AUPRC; seconds |
| `phecode-null` | logistic | `snr = 0` | AUROC ≈ 0.5 (calibration) |
| `phecode-oracle` | threshold | `γ = sd_pop = closure = parent_dropout = noise_frac = 0`, `λ = 1`, `snr = ∞`, dense retrieval, clean ICD layer | AUROC = AUPRC = 1 (calibration) |
| `phecode-xu` | log-linear | density from 5k patients × 100 events, κ = 1.5 (mean degree ~14) | KESER-like graph; ~1–3 min |

Default parameters were tuned to match a summary-only profile of real inputs
(mean degree, hierarchy compliance, k, score range)

## 5. Outputs

Each run writes to `<out>` (default `data_sim/<preset>/`):

| file | represents |
|---|---|
| `dim/Phecode_map_v1_2_icd9_icd10cm.csv`, `dim/phecode_definitions1.2.csv` | public Phecode 1.2 (real) |
| `dim/icd_athena_2025/CONCEPT.csv` | ICD typing dictionary (synthetic, from the public map) |
| `validation/dict-disease-disease-wiki-bin-eval.Rdata` | gold-standard dx–dx graph |
| `prediction/<query>_Diagnosis.csv` | ranked raw-ICD retrievals per query |
| `spec/{nodes,gold_edges,queries,predictions,embeddings}.csv` | phecode-level spec tables (`docs/SPEC-inputs.md`) |
| `sim_manifest.json`, `README.md` | every parameter, seed, source checksums (sha256s), per-query truth (`n_pos`, `recall_at_k`) |

## 6. Calibration tests these enable

- **Null**: `snr = 0` ⇒ AUROC ∈ [0.45, 0.55] and AUPRC ≈ prevalence in every
  cell (averaged over seeds).
- **Monotone**: AUROC and AUPRC non-decreasing in `snr` and in `λ`.
- **Oracle**: `λ = 1, snr = ∞, noise_frac = 0, k ≥ degree` ⇒ AUROC = AUPRC = 1;
  with `k < degree` the floored positives lower AUROC by exactly the expected
  amount.
- **Hierarchy-only baseline**: `λ = 0` ⇒ above chance in `random` and
  `non_adjacent`, near chance in `intracategory`, and the gap between them
  grows with `γ`.
- **Strata ordering** (soft): `non_adjacent ≥ random ≥ intracategory ≥
  degree_matched` on AUROC when `γ > 0`.
- **Sibling failure mode**: with `λ = 0`, false positives concentrate on
  sibling pairs absent from the gold graph; their share rises as `closure`
  falls.
- **Children-count bias** (via the ICD layer): `score_max` favours phecodes
  with many fine codes; `score_mean` does not.

## 7. Sampler invariant tests (independent of the generative model)

Exact cell sizes `min(n_pos × r, pool)`; nested cells across ratios; no gold
edge, self-loop, hierarchy pair, or duplicate among controls; category
constraints per stratum; query endpoint preserved in `degree_matched`;
determinism under `seed`; zero-positive query excluded; query–query positive
counted under both queries.
