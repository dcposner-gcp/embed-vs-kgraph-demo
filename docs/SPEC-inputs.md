# Spec tables — the input contract

The language-neutral, phecode-level view of each simulated dataset. Written
by `WriteSimSpecTables()` to `data_sim/<preset>/spec/` as five CSVs. The
pipeline itself consumes the raw `data/`-shaped tree (ICD-level retrievals,
`.Rdata` gold graph); these tables carry the same information
post-harmonisation, so an independent implementation that reads them should
reproduce the pipeline's metrics.

## 1.1 `nodes.csv` — the vocabulary

| column | type | contents |
|---|---|---|
| `node_id` | string | unique; the Phecode 1.2 code as written (e.g. `250.2`) |
| `label` | string | phecode description |
| `category` | string | phecode group name; drives the `non_adjacent` / `intracategory` control strata |
| `parent_id` | string | dotted parent when it exists in the vocabulary; empty for roots |
| `node_type` | string | constant `dx` (single-type task; not used by the pipeline) |

## 1.2 `gold_edges.csv` — the knowledge-graph gold standard

| column | type | contents |
|---|---|---|
| `u`, `v` | string | both exist in `nodes`; `u < v` (string order), no self-loops, no ancestor/descendant pairs; undirected |
| `source` | string | provenance tag: `base` (edge rule), `closure` (surviving hierarchy-inherited edge), `noise` (uniform rewiring), `topup` (query min-degree fill), `forced` (forced query–query pair) |

## 1.3 `queries.csv` — the query set

| column | type | contents |
|---|---|---|
| `query_id` | string | unique; the query's name (e.g. `T2DM`) |
| `label` | string | display name |
| `node_id` | string | the query's Phecode 1.2 code; **empty** for COVID and Long COVID (no canonical 1.2 code) — such queries are excluded from evaluation, though their raw retrievals still exist under `prediction/` |

## 1.4 `predictions.csv` — what the retriever returned

One row per (query, candidate), already in the gold vocabulary (the raw-ICD
form of the same retrievals, with lists/ranges/ambiguity, is what the
pipeline harmonises from `prediction/<query>_Diagnosis.csv`).

| column | type | contents |
|---|---|---|
| `query_id` | string | exists in `queries` |
| `candidate_id` | string | exists in `nodes`; the query's own self code is excluded |
| `rank` | integer | ≥ 1; unique within `query_id` |
| `score` | float | in (0, 1], higher = more related; a rank-preserving map of cosine, `(cos + 1) / 2` |

Retrieval is top-`k` with `k` drawn per query (`k_mean`, `k_sdlog`);
`k_mean = Inf` (the `phecode-oracle` preset) makes it dense. The pipeline's
floor policy — pairs a query never retrieved are treated as predicted
negatives — is what makes top-k and dense runs comparable.

## 1.5 `embeddings.csv` — the retriever's vectors

| column | type | contents |
|---|---|---|
| `node_id` | string | exists in `nodes` |
| `dim_1` … `dim_d` | float | the retriever's embedding `U` (`d = 32` by default), rounded to 6 dp |

Not consumed by the pipeline — provided so the prediction scores can be
reproduced (`score = (cos(U_q, U_j) + 1) / 2`, then top-`k`) or an
alternative scorer can be run against the identical gold standard and
controls.
