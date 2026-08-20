# Compatibility patches for running Monocle 2 with recent dplyr/igraph.






.patch_binding <- function(namespace, name, value) {
  if (!exists(name, envir = namespace, inherits = FALSE)) {
    return(FALSE)
  }

  if (bindingIsLocked(name, namespace)) {
    unlockBinding(name, namespace)
  }

  assign(name, value, envir = namespace)
  lockBinding(name, namespace)

  TRUE
}

.scoped_verb_vars <- function(dots) {
  vars <- unlist(lapply(dots, function(x) {
    if (is.symbol(x)) {
      return(as.character(x))
    }

    x <- as.character(x)
    x[length(x)]
  }), use.names = FALSE)

  vars[!is.na(vars) & nzchar(vars)]
}

.parse_scoped_expressions <- function(dots) {
  exprs <- unlist(lapply(dots, as.character), use.names = FALSE)
  exprs <- exprs[!is.na(exprs) & nzchar(exprs)]
  rlang::parse_exprs(exprs)
}

.monocle2_select_ <- function(.data, ..., .dots = list()) {
  dots <- c(list(...), .dots)

  if (length(dots) == 0) {
    return(.data)
  }

  old_names <- .scoped_verb_vars(dots)
  old_names <- old_names[old_names %in% colnames(.data)]

  if (length(old_names) == 0) {
    return(.data[, character(0), drop = FALSE])
  }

  new_names <- names(dots)
  if (is.null(new_names)) {
    new_names <- old_names
  } else {
    new_names <- new_names[seq_along(old_names)]
    new_names[is.na(new_names) | new_names == ""] <- old_names[is.na(new_names) | new_names == ""]
  }

  out <- .data[, old_names, drop = FALSE]
  colnames(out) <- new_names
  out
}

.monocle2_group_by_ <- function(.data, ..., .dots = list()) {
  dots <- c(list(...), .dots)
  vars <- .scoped_verb_vars(dots)
  vars <- vars[vars %in% colnames(.data)]

  if (length(vars) == 0) {
    return(dplyr::group_by(.data))
  }

  dplyr::group_by(.data, !!!rlang::syms(vars))
}

.monocle2_mutate_ <- function(.data, ..., .dots = list()) {
  dots <- c(list(...), .dots)
  exprs <- .parse_scoped_expressions(dots)
  dplyr::mutate(.data, !!!exprs)
}

.monocle2_summarise_ <- function(.data, ..., .dots = list()) {
  dots <- c(list(...), .dots)
  exprs <- .parse_scoped_expressions(dots)
  dplyr::summarise(.data, !!!exprs)
}

.monocle2_filter_ <- function(.data, ..., .dots = list()) {
  dots <- c(list(...), .dots)
  exprs <- .parse_scoped_expressions(dots)
  dplyr::filter(.data, !!!exprs)
}

patch_dplyr_for_monocle2 <- function() {
  dplyr_ns <- asNamespace("dplyr")
  monocle_ns <- asNamespace("monocle")

  scoped_verbs <- list(
    select_ = .monocle2_select_,
    group_by_ = .monocle2_group_by_,
    mutate_ = .monocle2_mutate_,
    summarise_ = .monocle2_summarise_,
    filter_ = .monocle2_filter_
  )

  for (name in names(scoped_verbs)) {
    .patch_binding(dplyr_ns, name, scoped_verbs[[name]])
    .patch_binding(monocle_ns, name, scoped_verbs[[name]])
  }

  invisible(TRUE)
}

patch_igraph_dfs_for_monocle2 <- function() {
  igraph_ns <- asNamespace("igraph")
  monocle_ns <- asNamespace("monocle")

  old_dfs <- get("dfs", envir = igraph_ns)
  if (isTRUE(attr(old_dfs, "monocle2_patch"))) {
    return(invisible(TRUE))
  }

  patched_dfs <- function(graph, root, mode = "out", ..., neimode = NULL) {
    if (!is.null(neimode)) {
      mode <- neimode
    }

    old_dfs(graph = graph, root = root, mode = mode, ...)
  }

  attr(patched_dfs, "monocle2_patch") <- TRUE

  .patch_binding(igraph_ns, "dfs", patched_dfs)
  .patch_binding(monocle_ns, "dfs", patched_dfs)

  invisible(TRUE)
}

patch_project2MST_monocle2 <- function() {
  monocle_ns <- asNamespace("monocle")

  if (!exists("project2MST", envir = monocle_ns, inherits = FALSE)) {
    return(invisible(FALSE))
  }

  project2MST <- get("project2MST", envir = monocle_ns)
  if (isTRUE(attr(project2MST, "monocle2_patch"))) {
    return(invisible(TRUE))
  }

  txt <- paste(deparse(body(project2MST)), collapse = "\n")

  nei_patterns <- c(
    "names\\(V\\(dp_mst\\)\\[suppressWarnings\\(nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\)\\]\\)",
    "names\\(V\\(dp_mst\\)\\[suppressWarnings\\(igraph::\\.nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\)\\]\\)",
    "names\\(V\\(dp_mst\\)\\[suppressWarnings\\(\\.nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\)\\]\\)",
    "names\\(V\\(dp_mst\\)\\[igraph::\\.nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\]\\)",
    "names\\(V\\(dp_mst\\)\\[\\.nei\\(closest_vertex_names\\[i\\],\\s*mode\\s*=\\s*\"all\"\\)\\]\\)"
  )

  for (pattern in nei_patterns) {
    txt <- gsub(
      pattern,
      "names(igraph::neighbors(dp_mst, closest_vertex_names[i], mode = \"all\"))",
      txt,
      perl = TRUE
    )
  }

  txt <- gsub(
    "graph\\.adjacency\\(",
    "igraph::graph_from_adjacency_matrix(",
    txt,
    fixed = FALSE
  )

  txt <- gsub(
    "minimum\\.spanning\\.tree\\(",
    "igraph::mst(",
    txt,
    fixed = FALSE
  )

  body(project2MST) <- parse(text = txt)[[1]]
  environment(project2MST) <- monocle_ns
  attr(project2MST, "monocle2_patch") <- TRUE

  .patch_binding(monocle_ns, "project2MST", project2MST)

  invisible(TRUE)
}

patch_plot_cell_trajectory_monocle2 <- function() {
  monocle_ns <- asNamespace("monocle")

  if (!exists("plot_cell_trajectory", envir = monocle_ns, inherits = FALSE)) {
    return(invisible(FALSE))
  }

  plot_cell_trajectory <- get("plot_cell_trajectory", envir = monocle_ns)
  if (isTRUE(attr(plot_cell_trajectory, "monocle2_patch"))) {
    return(invisible(TRUE))
  }

  txt <- paste(deparse(body(plot_cell_trajectory)), collapse = "\n")

  txt <- gsub(
    "select_\\(prin_graph_dim_1 = x,\\s*prin_graph_dim_2 = y\\)",
    "dplyr::select(prin_graph_dim_1 = dplyr::all_of(x), prin_graph_dim_2 = dplyr::all_of(y))",
    txt,
    perl = TRUE
  )

  txt <- gsub(
    "select_\\(source = \"from\",\\s*target = \"to\"\\)",
    "dplyr::select(source = dplyr::all_of(\"from\"), target = dplyr::all_of(\"to\"))",
    txt,
    perl = TRUE
  )

  txt <- gsub(
    "select_\\(source = \"sample_name\",\\s*source_prin_graph_dim_1 = \"prin_graph_dim_1\",\\s*source_prin_graph_dim_2 = \"prin_graph_dim_2\"\\)",
    "dplyr::select(source = dplyr::all_of(\"sample_name\"), source_prin_graph_dim_1 = dplyr::all_of(\"prin_graph_dim_1\"), source_prin_graph_dim_2 = dplyr::all_of(\"prin_graph_dim_2\"))",
    txt,
    perl = TRUE
  )

  txt <- gsub(
    "select_\\(target = \"sample_name\",\\s*target_prin_graph_dim_1 = \"prin_graph_dim_1\",\\s*target_prin_graph_dim_2 = \"prin_graph_dim_2\"\\)",
    "dplyr::select(target = dplyr::all_of(\"sample_name\"), target_prin_graph_dim_1 = dplyr::all_of(\"prin_graph_dim_1\"), target_prin_graph_dim_2 = dplyr::all_of(\"prin_graph_dim_2\"))",
    txt,
    perl = TRUE
  )

  txt <- gsub(
    "select_\\(data_dim_1 = x,\\s*data_dim_2 = y\\)",
    "dplyr::select(data_dim_1 = dplyr::all_of(x), data_dim_2 = dplyr::all_of(y))",
    txt,
    perl = TRUE
  )

  body(plot_cell_trajectory) <- parse(text = txt)[[1]]
  environment(plot_cell_trajectory) <- monocle_ns
  attr(plot_cell_trajectory, "monocle2_patch") <- TRUE

  .patch_binding(monocle_ns, "plot_cell_trajectory", plot_cell_trajectory)
  patch_dplyr_for_monocle2()

  invisible(TRUE)
}

patch_monocle2_igraph_compat <- function() {
  patch_igraph_dfs_for_monocle2()
  patch_project2MST_monocle2()
  invisible(TRUE)
}

patch_monocle2 <- function() {
  patch_dplyr_for_monocle2()
  patch_monocle2_igraph_compat()
  patch_plot_cell_trajectory_monocle2()
  invisible(TRUE)
}

# Backward-compatible aliases used in older notebooks/scripts.
patch_dplyr_for_monocle2_v2 <- patch_dplyr_for_monocle2
patch_dplyr_select_rename_monocle2 <- patch_dplyr_for_monocle2
patch_monocle_select_ <- patch_dplyr_for_monocle2
patch_plot_cell_trajectory_monocle2_v2 <- patch_plot_cell_trajectory_monocle2
patch_plot_cell_trajectory_all_select <- patch_plot_cell_trajectory_monocle2
patch_plot_cell_trajectory_remaining_select <- patch_plot_cell_trajectory_monocle2
patch_igraph_for_monocle2 <- patch_monocle2_igraph_compat
patch_igraph_nei_for_monocle2 <- patch_monocle2_igraph_compat
patch_monocle2_nei <- patch_monocle2_igraph_compat
patch_project2MST_monocle2_v2 <- patch_project2MST_monocle2
