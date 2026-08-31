# embed-vs-kgraph-demo — host-R wrappers. SIM_PRESET: phecode | phecode-null |
# phecode-oracle | phecode-xu; SIM_ARGS: extra simulate_data.R flags.
SIM_PRESET ?= phecode
SIM_ARGS   ?=

.PHONY: restore simulate pipeline-sim smoke

restore: ## one-time: restore the locked R library from renv.lock
	Rscript -e 'renv::restore(prompt = FALSE)'

simulate: ## generate simulated inputs -> data_sim/$(SIM_PRESET)/
	Rscript scripts/simulate_data.R --preset $(SIM_PRESET) --force $(SIM_ARGS)

pipeline-sim: ## run the unchanged pipeline -> scratch/sim_run/$(SIM_PRESET)/output/dx/
	./scripts/run_pipeline_sim.sh data_sim/$(SIM_PRESET) $(SIM_PRESET)

smoke: ## calibration checks: null -> AUROC ~ 0.5, oracle -> AUROC = 1
	$(MAKE) simulate pipeline-sim SIM_PRESET=phecode-null
	$(MAKE) simulate pipeline-sim SIM_PRESET=phecode-oracle
