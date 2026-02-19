# Master benchmark runner
# Sources all benchmark scripts, saves results as .rds, and optionally renders an HTML report.
# Usage: source("inst/benchmarks/run_all_benchmarks.R")

if (!requireNamespace("bench", quietly = TRUE))
  stop("Package 'bench' is required. Install with: install.packages('bench')")

bench_dir  <- file.path("inst", "benchmarks")
output_dir <- file.path(bench_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- Estimation benchmarks ----
message("Running estimation benchmarks...")
local({
  source(file.path(bench_dir, "bench_estimation.R"), local = TRUE)
  saveRDS(results_mle,   file.path(output_dir, paste0(timestamp, "_mle.rds")))
  saveRDS(results_ktee,  file.path(output_dir, paste0(timestamp, "_ktee.rds")))
})
message("  Saved to ", output_dir)

# ---- Simulation benchmarks ----
message("Running simulation benchmarks...")
local({
  source(file.path(bench_dir, "bench_simulation.R"), local = TRUE)
  saveRDS(results_sim,   file.path(output_dir, paste0(timestamp, "_sim.rds")))
  saveRDS(results_nsims, file.path(output_dir, paste0(timestamp, "_nsims.rds")))
  saveRDS(results_cond,  file.path(output_dir, paste0(timestamp, "_cond.rds")))
})
message("  Saved to ", output_dir)

# ---- C++ function benchmarks ----
message("Running C++ benchmarks...")
local({
  source(file.path(bench_dir, "bench_cpp.R"), local = TRUE)
  saveRDS(results_md,   file.path(output_dir, paste0(timestamp, "_md.rds")))
  saveRDS(results_hd,   file.path(output_dir, paste0(timestamp, "_hd.rds")))
  saveRDS(results_draw, file.path(output_dir, paste0(timestamp, "_draw.rds")))
  saveRDS(results_wtp,  file.path(output_dir, paste0(timestamp, "_wtp.rds")))
  saveRDS(results_cmd,  file.path(output_dir, paste0(timestamp, "_cmd.rds")))
})
message("  Saved to ", output_dir)

message("\nAll benchmarks complete. Results saved to: ", output_dir)
message("To load a result: readRDS('", output_dir, "/", timestamp, "_mle.rds')")

# ---- Optional HTML report ----
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  rmd_file <- file.path(bench_dir, "benchmark_report.Rmd")
  if (file.exists(rmd_file)) {
    report_out <- file.path(output_dir, paste0(timestamp, "_report.html"))
    message("\nRendering HTML report to: ", report_out)
    rmarkdown::render(rmd_file, output_file = report_out, quiet = TRUE)
    message("  Done.")
  }
}
