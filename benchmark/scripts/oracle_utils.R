# oracle_utils_codex_20260331.R
# Graph-only oracle utilities for Benchmark V1
# HARD RULE: this file must not use splits/bipartitions to classify missingness.

suppressPackageStartupMessages({
  library(ape)
})

ROOT_SENTINEL <- "__ROOT__"

freeze_full_tree_identity <- function(full_tree,
                                      tree_semantics = c("unrooted", "rooted")) {
  stopifnot(inherits(full_tree, "phylo"))
  tree_semantics <- match.arg(tree_semantics)

  n_tip <- length(full_tree$tip.label)
  n_node <- full_tree$Nnode

  full_tree$node.label <- paste0("N_", (n_tip + 1):(n_tip + n_node))

  edge_mat <- full_tree$edge
  edge_len <- full_tree$edge.length
  if (is.null(edge_len)) {
    edge_len <- rep(NA_real_, nrow(edge_mat))
  }

  get_node_label <- function(node_id) {
    if (node_id <= n_tip) {
      full_tree$tip.label[node_id]
    } else {
      full_tree$node.label[node_id - n_tip]
    }
  }

  parent_labels <- vapply(edge_mat[, 1], get_node_label, character(1))
  child_labels <- vapply(edge_mat[, 2], get_node_label, character(1))
  root_label <- setdiff(parent_labels, child_labels)[1]

  identity_tbl <- data.frame(
    edge_row = seq_len(nrow(edge_mat)),
    parent_node = edge_mat[, 1],
    child_node = edge_mat[, 2],
    parent_label = parent_labels,
    child_label = child_labels,
    branch_id = child_labels,
    branch_type = ifelse(edge_mat[, 2] <= n_tip, "terminal", "internal"),
    branch_length = edge_len,
    stringsAsFactors = FALSE
  )

  list(
    full_tree = full_tree,
    identity_tbl = identity_tbl,
    root_label = root_label,
    tree_semantics = tree_semantics
  )
}

build_descendant_branch_table <- function(full_tree) {
  stopifnot(inherits(full_tree, "phylo"))

  n_tip <- length(full_tree$tip.label)
  if (is.null(full_tree$node.label)) {
    full_tree$node.label <- paste0("N_", (n_tip + 1):(n_tip + full_tree$Nnode))
  }

  get_descendant_tips_by_node <- function(node_id) {
    if (node_id <= n_tip) {
      return(full_tree$tip.label[node_id])
    }
    kids <- full_tree$edge[full_tree$edge[, 1] == node_id, 2]
    sort(as.character(unlist(lapply(kids, get_descendant_tips_by_node), use.names = FALSE)))
  }

  edge_tbl <- freeze_full_tree_identity(full_tree)$identity_tbl
  out <- vector("list", nrow(edge_tbl))

  for (i in seq_len(nrow(edge_tbl))) {
    child_node <- edge_tbl$child_node[i]
    desc_tips <- get_descendant_tips_by_node(child_node)
    out[[i]] <- data.frame(
      branch_id = edge_tbl$branch_id[i],
      branch_type = edge_tbl$branch_type[i],
      stringsAsFactors = FALSE
    )
    out[[i]]$desc_tip_labels <- list(desc_tips)
  }

  tbl <- do.call(rbind, out)
  rownames(tbl) <- NULL
  tbl
}

get_plot_order <- function(tree) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)
  plot.phylo(tree, show.tip.label = TRUE, no.margin = TRUE)
  pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  yy <- pp$yy[seq_len(length(tree$tip.label))]
  names(yy) <- tree$tip.label
  names(sort(yy, decreasing = TRUE))
}

make_pruning_schedule <- function(tip_order_labels,
                                  schedule_type = c("global_ordered", "random", "local_forward", "local_backward"),
                                  start_rank = 1L,
                                  delete_count = NULL,
                                  delete_prop = NULL,
                                  min_retained_tip_count = 3L,
                                  max_delete_fraction = 0.7,
                                  random_seed = NULL) {
  schedule_type <- match.arg(schedule_type)
  tip_order_labels <- as.character(tip_order_labels)
  n <- length(tip_order_labels)

  if (!is.null(delete_count) && !is.null(delete_prop)) {
    stop("Use only one of delete_count or delete_prop")
  }
  if (n < min_retained_tip_count) {
    stop("tip count is smaller than min_retained_tip_count")
  }
  if (start_rank < 1L || start_rank > n) {
    stop("start_rank is out of range")
  }

  if (schedule_type == "random") {
    if (!is.null(random_seed)) {
      set.seed(random_seed)
    }
    ordered <- sample(tip_order_labels, size = n, replace = FALSE)
  } else {
    ordered <- tip_order_labels
  }

  max_by_fraction <- floor(n * max_delete_fraction)
  max_by_retained <- n - min_retained_tip_count

  direction_available_length <- switch(
    schedule_type,
    global_ordered = n,
    random = n,
    local_forward = n - start_rank + 1L,
    local_backward = start_rank
  )

  k_max <- min(direction_available_length, max_by_fraction, max_by_retained)

  if (is.null(delete_count) && is.null(delete_prop)) {
    delete_count <- k_max
  }
  if (!is.null(delete_prop)) {
    delete_count <- floor(n * delete_prop)
  }
  delete_count <- as.integer(max(0L, min(delete_count, k_max)))

  deleted_sequence <- switch(
    schedule_type,
    global_ordered = ordered[seq_len(delete_count)],
    random = ordered[seq_len(delete_count)],
    local_forward = ordered[start_rank:(start_rank + delete_count - 1L)],
    local_backward = rev(ordered[(start_rank - delete_count + 1L):start_rank])
  )
  if (delete_count == 0L) {
    deleted_sequence <- character(0)
  }

  out <- vector("list", delete_count + 1L)
  for (step_id in 0:delete_count) {
    deleted_tip_labels <- if (step_id == 0L) character(0) else deleted_sequence[seq_len(step_id)]
    deleted_tip_label <- if (step_id == 0L) NA_character_ else deleted_sequence[step_id]
    retained_tip_labels <- setdiff(tip_order_labels, deleted_tip_labels)

    out[[step_id + 1L]] <- data.frame(
      step_id = step_id,
      schedule_type = schedule_type,
      start_rank = start_rank,
      deleted_tip_label = deleted_tip_label,
      retained_tip_count = length(retained_tip_labels),
      max_delete_count = k_max,
      stringsAsFactors = FALSE
    )
    out[[step_id + 1L]]$deleted_tip_labels <- list(deleted_tip_labels)
    out[[step_id + 1L]]$retained_tip_labels <- list(retained_tip_labels)
  }

  schedule_tbl <- do.call(rbind, out)
  rownames(schedule_tbl) <- NULL
  schedule_tbl
}

build_subtree_series <- function(full_tree, schedule_tbl) {
  stopifnot(inherits(full_tree, "phylo"))

  subtree_list <- vector("list", nrow(schedule_tbl))
  names(subtree_list) <- paste0("step", schedule_tbl$step_id)

  for (i in seq_len(nrow(schedule_tbl))) {
    deleted <- schedule_tbl$deleted_tip_labels[[i]]
    subtree_list[[i]] <- if (length(deleted) == 0L) full_tree else drop.tip(full_tree, deleted)
  }

  list(subtree_list = subtree_list, subtree_meta_tbl = schedule_tbl)
}

init_graph_state <- function(identity_tbl, root_label,
                             tree_semantics = c("unrooted", "rooted")) {
  tree_semantics <- match.arg(tree_semantics)
  edges <- identity_tbl[, c("parent_label", "child_label", "branch_length", "branch_type")]
  edges$parent_label[edges$parent_label == root_label] <- ROOT_SENTINEL
  edges$member_ids <- lapply(identity_tbl$branch_id, function(x) x)
  edges$edge_uid <- seq_len(nrow(edges))
  rownames(edges) <- NULL

  list(
    edges = edges,
    next_edge_uid = nrow(edges) + 1L,
    tree_semantics = tree_semantics,
    tip_labels = as.character(identity_tbl$branch_id[identity_tbl$branch_type == "terminal"])
  )
}

member_string <- function(member_ids) {
  ids <- sort(unique(as.character(member_ids)))
  if (length(ids) == 0L) {
    return(NA_character_)
  }
  paste(ids, collapse = "|")
}

format_branch_length_value <- function(x) {
  if (is.na(x)) {
    return("")
  }
  format(x, digits = 10, scientific = FALSE, trim = TRUE)
}

find_child_edge <- function(state, child_label) {
  which(state$edges$child_label == child_label)
}

find_child_edges_of_parent <- function(state, parent_label) {
  which(state$edges$parent_label == parent_label)
}

find_incident_edges_of_node <- function(state, node_label) {
  which(state$edges$parent_label == node_label | state$edges$child_label == node_label)
}

other_endpoint <- function(state, row_idx, node_label) {
  if (state$edges$parent_label[row_idx] == node_label) {
    return(as.character(state$edges$child_label[row_idx]))
  }
  if (state$edges$child_label[row_idx] == node_label) {
    return(as.character(state$edges$parent_label[row_idx]))
  }
  stop(sprintf("Node %s is not incident to row %d", node_label, row_idx))
}

promote_root_singleton <- function(state, parent_label) {
  child_rows <- find_child_edges_of_parent(state, parent_label)
  if (length(child_rows) == 1L) {
    state$edges$parent_label[child_rows] <- ROOT_SENTINEL
  }
  state
}

merge_rows_into_edge <- function(state, row_a, row_b, parent_label, child_label) {
  len_a <- state$edges$branch_length[row_a]
  len_b <- state$edges$branch_length[row_b]
  merged_length <- if (all(is.na(c(len_a, len_b)))) NA_real_ else sum(c(len_a, len_b), na.rm = TRUE)
  merged_members <- sort(unique(c(
    state$edges$member_ids[[row_a]],
    state$edges$member_ids[[row_b]]
  )))
  child_type <- if (state$edges$branch_type[row_a] == "internal" ||
                    state$edges$branch_type[row_b] == "internal") {
    "internal"
  } else {
    "terminal"
  }

  drop_rows <- sort(c(row_a, row_b), decreasing = TRUE)
  for (r in drop_rows) {
    state$edges <- state$edges[-r, , drop = FALSE]
  }

  state$edges <- rbind(
    state$edges,
    data.frame(
      parent_label = parent_label,
      child_label = child_label,
      branch_length = merged_length,
      branch_type = child_type,
      member_ids = I(list(merged_members)),
      edge_uid = state$next_edge_uid,
      stringsAsFactors = FALSE
    )
  )
  state$next_edge_uid <- state$next_edge_uid + 1L
  rownames(state$edges) <- NULL
  state
}

choose_edge_orientation <- function(state, endpoint_a, endpoint_b) {
  endpoints <- c(as.character(endpoint_a), as.character(endpoint_b))

  if (ROOT_SENTINEL %in% endpoints) {
    other <- setdiff(endpoints, ROOT_SENTINEL)
    return(list(parent_label = ROOT_SENTINEL, child_label = other[1]))
  }

  is_tip <- endpoints %in% state$tip_labels
  if (sum(is_tip) == 1L) {
    return(list(
      parent_label = endpoints[!is_tip][1],
      child_label = endpoints[is_tip][1]
    ))
  }

  endpoints <- sort(endpoints)
  list(parent_label = endpoints[1], child_label = endpoints[2])
}

normalize_root_unrooted <- function(state) {
  if (!identical(state$tree_semantics, "unrooted")) {
    return(state)
  }

  repeat {
    root_rows <- find_child_edges_of_parent(state, ROOT_SENTINEL)
    if (length(root_rows) <= 1L) {
      state <- promote_root_singleton(state, ROOT_SENTINEL)
      break
    }
    if (length(root_rows) != 2L) {
      break
    }

    root_edges <- state$edges[root_rows, , drop = FALSE]
    ord <- order(root_edges$child_label, root_edges$edge_uid)
    root_rows <- root_rows[ord]

    parent_side <- root_rows[1]
    child_side <- root_rows[2]
    parent_label <- state$edges$child_label[parent_side]
    child_label <- state$edges$child_label[child_side]

    state <- merge_rows_into_edge(
      state = state,
      row_a = parent_side,
      row_b = child_side,
      parent_label = parent_label,
      child_label = child_label
    )
  }

  state
}

contract_degree2_nodes_unrooted <- function(state) {
  if (!identical(state$tree_semantics, "unrooted")) {
    return(state)
  }

  repeat {
    node_labels <- sort(unique(c(state$edges$parent_label, state$edges$child_label)))
    node_labels <- setdiff(node_labels, c(ROOT_SENTINEL, state$tip_labels))

    contracted <- FALSE
    for (node_label in node_labels) {
      incident_rows <- find_incident_edges_of_node(state, node_label)
      if (length(incident_rows) != 2L) {
        next
      }

      endpoint_a <- other_endpoint(state, incident_rows[1], node_label)
      endpoint_b <- other_endpoint(state, incident_rows[2], node_label)
      orient <- choose_edge_orientation(state, endpoint_a, endpoint_b)

      state <- merge_rows_into_edge(
        state = state,
        row_a = incident_rows[1],
        row_b = incident_rows[2],
        parent_label = orient$parent_label,
        child_label = orient$child_label
      )
      state <- normalize_root_unrooted(state)
      contracted <- TRUE
      break
    }

    if (!contracted) {
      break
    }
  }

  state
}

delete_tip_once_rooted <- function(state, deleted_tip_label) {
  tip_row <- find_child_edge(state, deleted_tip_label)
  if (length(tip_row) != 1L) {
    stop(sprintf("Expected exactly one current tip edge for %s", deleted_tip_label))
  }

  parent_label <- state$edges$parent_label[tip_row]
  sibling_rows <- setdiff(find_child_edges_of_parent(state, parent_label), tip_row)
  up_row <- find_child_edge(state, parent_label)

  local_e_term <- deleted_tip_label
  local_e_sib <- if (length(sibling_rows) == 1L) member_string(state$edges$member_ids[[sibling_rows]]) else NA_character_
  local_e_up <- if (length(up_row) == 1L) member_string(state$edges$member_ids[[up_row]]) else NA_character_
  local_sib_child_label <- if (length(sibling_rows) == 1L) as.character(state$edges$child_label[sibling_rows]) else NA_character_
  local_sib_parent_label <- if (length(sibling_rows) == 1L) as.character(state$edges$parent_label[sibling_rows]) else NA_character_
  local_up_child_label <- if (length(up_row) == 1L) as.character(state$edges$child_label[up_row]) else NA_character_
  local_up_parent_label <- if (length(up_row) == 1L) as.character(state$edges$parent_label[up_row]) else NA_character_

  state$edges <- state$edges[-tip_row, , drop = FALSE]

  sibling_rows <- find_child_edges_of_parent(state, parent_label)
  up_row <- find_child_edge(state, parent_label)

  if (length(sibling_rows) == 1L && length(up_row) == 1L) {
    sib_row <- sibling_rows[1]
    parent_of_parent <- state$edges$parent_label[up_row]
    child_label <- state$edges$child_label[sib_row]
    state <- merge_rows_into_edge(
      state = state,
      row_a = sib_row,
      row_b = up_row,
      parent_label = parent_of_parent,
      child_label = child_label
    )
  } else if (length(sibling_rows) == 1L && identical(parent_label, ROOT_SENTINEL)) {
    state <- promote_root_singleton(state, ROOT_SENTINEL)
  }

  rownames(state$edges) <- NULL
  list(
    state = state,
    event = list(
      deleted_tip_parent_label = as.character(parent_label),
      local_e_term = local_e_term,
      local_e_sib = local_e_sib,
      local_e_up = local_e_up,
      local_sib_child_label = local_sib_child_label,
      local_sib_parent_label = local_sib_parent_label,
      local_up_child_label = local_up_child_label,
      local_up_parent_label = local_up_parent_label
    )
  )
}

delete_tip_once_unrooted <- function(state, deleted_tip_label) {
  tip_rows <- find_incident_edges_of_node(state, deleted_tip_label)
  if (length(tip_rows) != 1L) {
    stop(sprintf("Expected exactly one current incident edge for tip %s", deleted_tip_label))
  }

  tip_row <- tip_rows[1]
  parent_label <- other_endpoint(state, tip_row, deleted_tip_label)
  incident_rows <- setdiff(find_incident_edges_of_node(state, parent_label), tip_row)

  local_e_term <- deleted_tip_label
  local_e_sib <- if (length(incident_rows) >= 1L) member_string(state$edges$member_ids[[incident_rows[1]]]) else NA_character_
  local_e_up <- if (length(incident_rows) >= 2L) member_string(state$edges$member_ids[[incident_rows[2]]]) else NA_character_
  local_sib_child_label <- if (length(incident_rows) >= 1L) as.character(state$edges$child_label[incident_rows[1]]) else NA_character_
  local_sib_parent_label <- if (length(incident_rows) >= 1L) as.character(state$edges$parent_label[incident_rows[1]]) else NA_character_
  local_up_child_label <- if (length(incident_rows) >= 2L) as.character(state$edges$child_label[incident_rows[2]]) else NA_character_
  local_up_parent_label <- if (length(incident_rows) >= 2L) as.character(state$edges$parent_label[incident_rows[2]]) else NA_character_

  state$edges <- state$edges[-tip_row, , drop = FALSE]
  state <- normalize_root_unrooted(state)
  state <- contract_degree2_nodes_unrooted(state)

  rownames(state$edges) <- NULL
  list(
    state = state,
    event = list(
      deleted_tip_parent_label = as.character(parent_label),
      local_e_term = local_e_term,
      local_e_sib = local_e_sib,
      local_e_up = local_e_up,
      local_sib_child_label = local_sib_child_label,
      local_sib_parent_label = local_sib_parent_label,
      local_up_child_label = local_up_child_label,
      local_up_parent_label = local_up_parent_label
    )
  )
}

delete_tip_once <- function(state, deleted_tip_label) {
  if (identical(state$tree_semantics, "unrooted")) {
    return(delete_tip_once_unrooted(state, deleted_tip_label))
  }
  delete_tip_once_rooted(state, deleted_tip_label)
}

state_to_classification <- function(state, identity_tbl) {
  axis_ids <- as.character(identity_tbl$branch_id)
  current_members <- state$edges$member_ids

  line <- rep(NA_character_, length(axis_ids))
  names(line) <- axis_ids

  merge_groups <- list()
  for (i in seq_len(nrow(state$edges))) {
    members <- sort(unique(as.character(current_members[[i]])))
    if (length(members) == 1L) {
      len_i <- state$edges$branch_length[i]
      line[members] <- format_branch_length_value(len_i)
    } else if (length(members) >= 2L) {
      line[members] <- "NA_fuse"
      merge_groups[[length(merge_groups) + 1L]] <- members
    }
  }

  still_na <- names(line)[is.na(line)]
  line[still_na] <- "NA_struct"

  list(
    line = line,
    merge_groups = merge_groups,
    na_struct = names(line)[line == "NA_struct"],
    na_fuse = names(line)[line == "NA_fuse"]
  )
}

classification_to_long_rows <- function(classification, run_id, step_id, semantics_label) {
  branch_ids <- names(classification$line)
  values <- unname(classification$line)
  statuses <- ifelse(
    values == "NA_struct",
    "NA_struct",
    ifelse(values == "NA_fuse", "NA_fuse", "observed")
  )

  data.frame(
    gene_id = paste0(run_id, "_step", step_id),
    step_id = step_id,
    branch_id = branch_ids,
    value = values,
    status = statuses,
    tree_semantics = semantics_label,
    stringsAsFactors = FALSE
  )
}

state_to_fusion_group_rows <- function(state,
                                       run_id,
                                       step_id,
                                       member_col_name = "benchmark_members") {
  rows <- list()
  group_counter <- 0L

  for (i in seq_len(nrow(state$edges))) {
    members <- sort(unique(as.character(state$edges$member_ids[[i]])))
    if (length(members) < 2L) {
      next
    }

    group_counter <- group_counter + 1L
    len_i <- state$edges$branch_length[i]
    row <- data.frame(
      gene_id = paste0(run_id, "_step", step_id),
      step_id = step_id,
      merge_group_id = sprintf("MG%02d", group_counter),
      expected_fused_length = format_branch_length_value(len_i),
      group_size = length(members),
      stringsAsFactors = FALSE
    )
    row[[member_col_name]] <- paste(members, collapse = "|")
    rows[[length(rows) + 1L]] <- row
  }

  if (length(rows) == 0L) {
    out <- data.frame(
      gene_id = character(),
      step_id = integer(),
      merge_group_id = character(),
      expected_fused_length = character(),
      group_size = integer(),
      stringsAsFactors = FALSE
    )
    out[[member_col_name]] <- character()
    return(out)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

active_edges_from_state <- function(state) {
  unique(unlist(state$edges$member_ids, use.names = FALSE))
}

trajectory_summary_from_classification <- function(classification, state) {
  sizes <- if (length(classification$merge_groups) == 0L) integer() else lengths(classification$merge_groups)
  data.frame(
    na_struct_count = length(classification$na_struct),
    na_fuse_count = length(classification$na_fuse),
    fusion_burden = if (length(sizes) == 0L) 0L else sum(sizes - 1L),
    max_group_size = if (length(sizes) == 0L) 1L else max(sizes),
    num_groups_ge2 = sum(sizes >= 2L),
    active_edge_count = length(active_edges_from_state(state)),
    stringsAsFactors = FALSE
  )
}
