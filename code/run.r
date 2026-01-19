
# -------------------------------
# Restore environment
# -------------------------------
message("Restoring environment with renv::restore()...")
getwd()
.libPaths()

if (!require(renv)) install.packages("renv", repos="https://cloud.r-project.org")

restore_success <- tryCatch(
  {
    renv::restore()
    TRUE
  },
  error = function(e) {
    message("renv::restore() failed: ", e$message)
    FALSE
  }
)


if (!restore_success) {
  message("Environment restoration failed. Exiting.")
  quit(status = 1)
}

library(ewaff)
library(lmerTest)

message("Environment successfully built. Data-dependent analysis not executed due to data restrictions.")