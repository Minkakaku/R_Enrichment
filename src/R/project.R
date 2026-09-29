# 合并公共配置和项目参数，显式切换到项目根目录后调用公共入口。
load_method <- function(project_dir, method) {
  project <- normalizePath(project_dir, mustWork = TRUE)
  workspace <- normalizePath(file.path(project, "../.."), mustWork = TRUE)
  config <- list()
  for (path in c(file.path(workspace, "config/default.yaml"),
                 file.path(workspace, "config/reference.yaml"),
                 file.path(project, "config/config.yaml"))) {
    config <- utils::modifyList(config, yaml::read_yaml(path), keep.null = TRUE)
  }
  options <- config$methods[[method]]
  path_keys <- c("input", "rds", "output", "object-dir", "object_dir", "gene-order", "database-dir", "gmt", "gsea-gmt")
  for (key in intersect(names(options), path_keys)) {
    value <- options[[key]]
    if (!is.null(value)) options[[key]] <- if (grepl("^/", path.expand(value))) path.expand(value) else file.path(project, value)
  }
  args <- unlist(lapply(names(options), function(key) {
    value <- options[[key]]
    if (is.null(value)) return(character())
    c(paste0("--", key), as.character(unlist(value)))
  }), use.names = FALSE)
  list(workspace = workspace, project = project, args = args)
}

run_project_method <- function(project_dir, method, extra = commandArgs(trailingOnly = TRUE)) {
  context <- load_method(project_dir, method)
  script <- file.path(context$workspace, "src/R", paste0(method, "_workflow.R"))
  args <- c(context$args, extra[extra != "--dry-run"])
  if ("--dry-run" %in% extra) {
    cat("cwd:", context$project, "\nscript:", script, "\nargs:", paste(shQuote(args), collapse = " "), "\n")
    return(invisible(context))
  }
  previous <- setwd(context$project)
  on.exit(setwd(previous), add = TRUE)
  dir.create("logs", showWarnings = FALSE)
  log <- file(file.path("logs", paste0(method, "_", format(Sys.time(), "%Y%m%d_%H%M%OS6"), ".log")), "wt")
  sink(log); sink(log, type = "message")
  on.exit({ sink(type = "message"); sink(); close(log) }, add = TRUE)
  cat("Arguments:", paste(args, collapse = " "), "\n")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)
  env$main(args)
}
