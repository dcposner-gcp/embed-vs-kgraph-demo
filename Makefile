# embed-vs-kgraph-demo — host-R wrappers. SIM_PRESET: phecode | phecode-null |
# phecode-oracle | phecode-xu; SIM_ARGS: extra simulate_data.R flags.
# dag / figures need GraphViz `dot` (figures also rsvg-convert for the PNG).
SIM_PRESET ?= phecode
SIM_ARGS   ?=

.PHONY: restore simulate pipeline-sim smoke dag figures

restore: ## one-time: restore the locked R library from renv.lock
	Rscript -e 'renv::restore(prompt = FALSE)'

simulate: ## generate simulated inputs -> data_sim/$(SIM_PRESET)/
	Rscript scripts/simulate_data.R --preset $(SIM_PRESET) --force $(SIM_ARGS)

pipeline-sim: ## run the unchanged pipeline -> scratch/sim_run/$(SIM_PRESET)/output/dx/
	./scripts/run_pipeline_sim.sh data_sim/$(SIM_PRESET) $(SIM_PRESET)

smoke: ## calibration checks: null -> AUROC ~ 0.5, oracle -> AUROC = 1
	$(MAKE) simulate pipeline-sim SIM_PRESET=phecode-null
	$(MAKE) simulate pipeline-sim SIM_PRESET=phecode-oracle

dag: ## render the targets DAG from _targets.R -> docs/figures/pipeline_dag.svg
	Rscript scripts/render_dag.R

figures: dag ## dag + the concept figure -> docs/figures/concept.{svg,png}
	docs/figures/concept.sh
