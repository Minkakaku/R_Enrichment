source("src/R/project.R")
project <- normalizePath("projects/project1")
method <- names(yaml::read_yaml("projects/project1/config/config.yaml")$methods)[1]
context <- load_method(project, method)
stopifnot(context$project == project, length(context$args) > 0)
for (script in list.files("src/R", pattern = "_workflow[.]R$", full.names = TRUE)) {
  env <- new.env()
  sys.source(script, envir = env)
  stopifnot(is.function(env$main))
}
run_project_method(project, method, extra = "--dry-run")
cat("Safe imports and method configuration checks passed\n")
