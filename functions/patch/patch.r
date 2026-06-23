patch_dplyr_for_monocle2 <- function() {
  ns <- asNamespace("dplyr")

  patch_fun <- function(name, fun) {
    if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
    assign(name, fun, envir = ns)
    lockBinding(name, ns)
  }

  patch_fun("select_", function(.data, ..., .dots = list()) {
    dots <- c(list(...), .dots)
    vars <- unlist(lapply(dots, as.character))
    dplyr::select(.data, dplyr::all_of(vars))
  })

  patch_fun("group_by_", function(.data, ..., .dots = list()) {
    dots <- c(list(...), .dots)
    vars <- unlist(lapply(dots, as.character))
    dplyr::group_by(.data, !!!rlang::syms(vars))
  })

  invisible(TRUE)
}


patch_igraph_for_monocle2 <- function() {
  ns <- asNamespace("igraph")

  if (bindingIsLocked("dfs", ns)) unlockBinding("dfs", ns)
  old_dfs <- get("dfs", envir = ns)

  assign("dfs", function(graph, root, mode = "out", ..., neimode = NULL) {
    if (!is.null(neimode)) mode <- neimode
    old_dfs(graph = graph, root = root, mode = mode, ...)
  }, envir = ns)

  lockBinding("dfs", ns)

  if (exists("nei", envir = ns, inherits = FALSE)) {
    if (bindingIsLocked("nei", ns)) unlockBinding("nei", ns)

    assign("nei", function(v, mode = "all") {
      igraph::.nei(v, mode = mode)
    }, envir = ns)

    lockBinding("nei", ns)
  }

  invisible(TRUE)
}


patch_igraph_nei_for_monocle2 <- function() {
  ns <- asNamespace("igraph")

  if (exists("nei", envir = ns, inherits = FALSE)) {
    if (bindingIsLocked("nei", ns)) unlockBinding("nei", ns)

    assign("nei", function(..., mode = "all") {
      igraph::.nei(..., mode = mode)
    }, envir = ns)

    lockBinding("nei", ns)
  }

  invisible(TRUE)
}


patch_monocle2_nei <- function() {
  ns <- asNamespace("monocle")

  if (exists("nei", envir = ns, inherits = FALSE)) {
    if (bindingIsLocked("nei", ns)) unlockBinding("nei", ns)

    assign("nei", function(v, mode = "all") {
      igraph::.nei(v, mode = mode)
    }, envir = ns)

    lockBinding("nei", ns)
  }

  invisible(TRUE)
}

patch_project2MST_monocle2 <- function() {
  ns <- asNamespace("monocle")

  f <- get("project2MST", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    "nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)",
    "igraph::.nei(closest_vertex_names[i], mode = \"all\")",
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("project2MST", ns)
  assign("project2MST", f, envir = ns)
  lockBinding("project2MST", ns)

  invisible(TRUE)
}


patch_project2MST_monocle2_v2 <- function() {
  ns <- asNamespace("monocle")

  f <- get("project2MST", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    "suppressWarnings\\(igraph::.nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\)",
    ".nei(closest_vertex_names[i], mode = \"all\")",
    txt
  )

  txt <- gsub(
    "suppressWarnings\\(nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\)",
    ".nei(closest_vertex_names[i], mode = \"all\")",
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("project2MST", ns)
  assign("project2MST", f, envir = ns)
  lockBinding("project2MST", ns)

  invisible(TRUE)
}

patch_dplyr_for_monocle2_v2 <- function() {
  ns <- asNamespace("dplyr")

  patch_fun <- function(name, fun) {
    if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
    assign(name, fun, envir = ns)
    lockBinding(name, ns)
  }

  patch_fun("select_", function(.data, ..., .dots = list()) {
    dots <- c(list(...), .dots)
    vars <- unlist(lapply(dots, function(x) {
      x <- as.character(x)
      x[length(x)]
    }))

    vars <- vars[vars %in% colnames(.data)]

    dplyr::select(.data, dplyr::all_of(vars))
  })

  patch_fun("group_by_", function(.data, ..., .dots = list()) {
    dots <- c(list(...), .dots)
    vars <- unlist(lapply(dots, function(x) {
      x <- as.character(x)
      x[length(x)]
    }))

    vars <- vars[vars %in% colnames(.data)]

    dplyr::group_by(.data, !!!rlang::syms(vars))
  })

  patch_fun("mutate_", function(.data, ..., .dots = list()) {
    dplyr::mutate(.data, ...)
  })

  patch_fun("summarise_", function(.data, ..., .dots = list()) {
    dplyr::summarise(.data, ...)
  })

  patch_fun("filter_", function(.data, ..., .dots = list()) {
    dplyr::filter(.data, ...)
  })

  invisible(TRUE)
}

patch_dplyr_select_rename_monocle2 <- function() {
  ns <- asNamespace("dplyr")

  if (bindingIsLocked("select_", ns)) unlockBinding("select_", ns)

  assign("select_", function(.data, ..., .dots = list()) {
    dots <- c(list(...), .dots)

    if (length(dots) == 0) {
      return(.data)
    }

    old_names <- unlist(lapply(dots, function(x) {
      as.character(x)[length(as.character(x))]
    }))

    new_names <- names(dots)

    if (is.null(new_names)) {
      new_names <- old_names
    }

    new_names[new_names == ""] <- old_names[new_names == ""]

    out <- dplyr::select(.data, dplyr::all_of(old_names))
    colnames(out) <- new_names

    out
  }, envir = ns)

  lockBinding("select_", ns)

  invisible(TRUE)
}


patch_monocle_select_ <- function() {
  ns <- asNamespace("monocle")

  if (exists("select_", envir = ns, inherits = FALSE)) {
    if (bindingIsLocked("select_", ns)) unlockBinding("select_", ns)

    assign("select_", function(.data, ..., .dots = list()) {
      dots <- c(list(...), .dots)

      old_names <- unlist(lapply(dots, function(x) {
        as.character(x)[length(as.character(x))]
      }))

      new_names <- names(dots)
      if (is.null(new_names)) new_names <- old_names
      new_names[new_names == ""] <- old_names[new_names == ""]

      out <- dplyr::select(.data, dplyr::all_of(old_names))
      colnames(out) <- new_names

      out
    }, envir = ns)

    lockBinding("select_", ns)
  }

  invisible(TRUE)
}


patch_plot_cell_trajectory_monocle2 <- function() {
  ns <- asNamespace("monocle")

  f <- get("plot_cell_trajectory", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    "select_\\(prin_graph_dim_1 = x,\\s*prin_graph_dim_2 = y\\)",
    "dplyr::select(., prin_graph_dim_1 = dplyr::all_of(x), prin_graph_dim_2 = dplyr::all_of(y))",
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("plot_cell_trajectory", ns)
  assign("plot_cell_trajectory", f, envir = ns)
  lockBinding("plot_cell_trajectory", ns)

  invisible(TRUE)
}

patch_plot_cell_trajectory_monocle2_v2 <- function() {
  ns <- asNamespace("monocle")

  f <- get("plot_cell_trajectory", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    "select_\\(prin_graph_dim_1 = x,\\s*\\n\\s*prin_graph_dim_2 = y\\)",
    "dplyr::select(., prin_graph_dim_1 = dplyr::all_of(x), prin_graph_dim_2 = dplyr::all_of(y))",
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("plot_cell_trajectory", ns)
  assign("plot_cell_trajectory", f, envir = ns)
  lockBinding("plot_cell_trajectory", ns)

  invisible(TRUE)
}


patch_plot_cell_trajectory_all_select <- function() {
  ns <- asNamespace("monocle")
  f <- get("plot_cell_trajectory", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    'select_\\(source = "from",\\s*target = "to"\\)',
    'dplyr::select(., source = "from", target = "to")',
    txt
  )

  txt <- gsub(
    'select_\\(source = "sample_name",\\s*prin_graph_dim_1,\\s*prin_graph_dim_2\\)',
    'dplyr::select(., source = "sample_name", prin_graph_dim_1, prin_graph_dim_2)',
    txt
  )

  txt <- gsub(
    'select_\\(target = "sample_name",\\s*prin_graph_dim_1,\\s*prin_graph_dim_2\\)',
    'dplyr::select(., target = "sample_name", prin_graph_dim_1, prin_graph_dim_2)',
    txt
  )

  txt <- gsub(
    'select_\\(data_dim_1 = x,\\s*data_dim_2 = y\\)',
    'dplyr::select(., data_dim_1 = dplyr::all_of(x), data_dim_2 = dplyr::all_of(y))',
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("plot_cell_trajectory", ns)
  assign("plot_cell_trajectory", f, envir = ns)
  lockBinding("plot_cell_trajectory", ns)

  invisible(TRUE)
}

patch_plot_cell_trajectory_remaining_select <- function() {
  ns <- asNamespace("monocle")
  f <- get("plot_cell_trajectory", envir = ns)

  txt <- paste(deparse(body(f)), collapse = "\n")

  txt <- gsub(
    'select_\\(source = "sample_name",\\s*source_prin_graph_dim_1 = "prin_graph_dim_1",\\s*source_prin_graph_dim_2 = "prin_graph_dim_2"\\)',
    'dplyr::select(., source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2")',
    txt
  )

  txt <- gsub(
    'select_\\(target = "sample_name",\\s*target_prin_graph_dim_1 = "prin_graph_dim_1",\\s*target_prin_graph_dim_2 = "prin_graph_dim_2"\\)',
    'dplyr::select(., target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2")',
    txt
  )

  body(f) <- parse(text = txt)[[1]]
  environment(f) <- ns

  unlockBinding("plot_cell_trajectory", ns)
  assign("plot_cell_trajectory", f, envir = ns)
  lockBinding("plot_cell_trajectory", ns)

  invisible(TRUE)
}
