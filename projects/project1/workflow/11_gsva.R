# 定位项目后调用共享方法；支持 --dry-run 检查配置。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project.R"))
run_project_method(project, "gsva")
