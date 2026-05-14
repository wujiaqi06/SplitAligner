#!/usr/bin/env Rscript

# benchmark.R
# Clean graph-only benchmark driver for Benchmark V1

suppressPackageStartupMessages({
  library(ape)
})

find_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  for (frm in rev(sys.frames())) {
    if (!is.null(frm$ofile)) {
      return(dirname(normalizePath(frm$ofile)))
    }
  }

  if (interactive()) {
    return(getwd())
  }

  stop("Unable to determine script directory")
}

script_dir <- find_script_dir()
base_dir <- dirname(script_dir)
find_oracle_utils <- function(dir_path) {
  candidates <- c(
    file.path(dir_path, "oracle_utils.R"),
    file.path(dir_path, "oracle_utils_codex_20260331.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop(sprintf("Cannot find oracle utils in %s", dir_path))
  }
  hit[1]
}
source(find_oracle_utils(script_dir))

get_env_or_default <- function(name, default = NULL) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    return(default)
  }
  val
}

get_env_int_or_default <- function(name, default = NULL) {
  val <- get_env_or_default(name, NULL)
  if (is.null(val)) {
    return(default)
  }
  as.integer(val)
}

get_env_num_or_default <- function(name, default = NULL) {
  val <- get_env_or_default(name, NULL)
  if (is.null(val)) {
    return(default)
  }
  as.numeric(val)
}

# -----------------------------
# User parameters (edit here)
# -----------------------------

SEED <- 42
N_TIPS_IF_SIM <- 10

# set TREE_FILE to a Newick path; otherwise a random tree is simulated
TREE_FILE <- get_env_or_default("TREE_FILE", file.path(base_dir, "inputs", "example_tree.nwk"))

# tree semantics:
#   unrooted / rooted
# Benchmark V1 defaults to unrooted semantics to match SplitAligner.
TREE_SEMANTICS <- get_env_or_default("TREE_SEMANTICS", "unrooted")

# packaged scenario name:
#   t10_global_deletion / t8_to_t3_local_deletion
SCENARIO_NAME <- get_env_or_default("SCENARIO_NAME", "t10_global_deletion")

# schedule types:
#   global_ordered / local_forward / local_backward / random
SCHEDULE_TYPE <- get_env_or_default("SCHEDULE_TYPE", "global_ordered")
START_RANK <- get_env_int_or_default("START_RANK", 1L)
DELETE_COUNT <- get_env_int_or_default("DELETE_COUNT", NULL)
DELETE_PROP <- get_env_num_or_default("DELETE_PROP", NULL)
MIN_RETAINED_TIPS <- get_env_int_or_default("MIN_RETAINED_TIPS", 3L)
# Benchmark V1 packaged default.
# This packaged demo uses 0.7 to make collapse / fusion progression more visible.
MAX_DELETE_FRACTION <- get_env_num_or_default("MAX_DELETE_FRACTION", 0.7)

# if > 0 and SCHEDULE_TYPE == "random", multiple replicates are run
NREP_RANDOM <- get_env_int_or_default("NREP_RANDOM", 0L)

SCENARIO_DIR <- file.path(base_dir, "outputs", SCENARIO_NAME)
SEMANTICS_SUBDIR <- if (identical(TREE_SEMANTICS, "rooted")) {
  "benchmark_rooted"
} else {
  "benchmark_unrooted"
}
OUT_DIR <- file.path(SCENARIO_DIR, SEMANTICS_SUBDIR)
OUT_PREFIX <- file.path(OUT_DIR, "benchmark")

scenario_presets <- list(
  t10_global_deletion = list(
    schedule_name = "fixed_t10_global_deletion",
    start_rank = 1L,
    deleted_sequence = c("t10", "t1", "t8", "t7", "t4", "t9", "t5")
  ),
  t8_to_t3_local_deletion = list(
    schedule_name = "fixed_middelete_t8_to_t3",
    start_rank = 3L,
    deleted_sequence = c("t8", "t7", "t4", "t9", "t5", "t2", "t3")
  )
)

write_species_tree_file <- function(tree, outfile) {
  writeLines(write.tree(tree), con = outfile)
}

write_gene_tree_file <- function(subtree_list, run_id, outfile) {
  con <- file(outfile, open = "wt")
  on.exit(close(con), add = TRUE)

  for (nm in names(subtree_list)) {
    step_id <- sub("^step", "", nm)
    gene_id <- paste0(run_id, "_", nm)
    tree_txt <- write.tree(subtree_list[[nm]])
    writeLines(paste0(gene_id, tree_txt), con = con)
  }
}

prepare_benchmark_tree <- function(tr,
                                   tree_semantics = c("unrooted", "rooted")) {
  tree_semantics <- match.arg(tree_semantics)
  if (is.null(tr$edge.length)) {
    stop("Benchmark tree must have branch lengths")
  }
  if (anyNA(tr$edge.length)) {
    stop("Benchmark tree contains missing branch lengths")
  }
  if (identical(tree_semantics, "unrooted") && is.rooted(tr)) {
    tr <- unroot(tr)
  }
  tr
}

read_or_simulate_tree <- function(tree_file, n_tips, seed, tree_semantics) {
  set.seed(seed)
  if (!is.null(tree_file) && file.exists(tree_file)) {
    tr <- read.tree(tree_file)
    if (inherits(tr, "multiPhylo")) stop("TREE_FILE must contain exactly one tree")
    return(prepare_benchmark_tree(tr, tree_semantics = tree_semantics))
  }
  prepare_benchmark_tree(rtree(n_tips), tree_semantics = tree_semantics)
}

make_fixed_pruning_schedule <- function(tip_order_labels,
                                        deleted_sequence,
                                        start_rank = 1L,
                                        schedule_name = "fixed_schedule",
                                        min_retained_tip_count = 3L) {
  tip_order_labels <- as.character(tip_order_labels)
  deleted_sequence <- as.character(deleted_sequence)

  if (anyDuplicated(deleted_sequence)) {
    stop("deleted_sequence contains duplicated tips")
  }
  if (!all(deleted_sequence %in% tip_order_labels)) {
    bad <- setdiff(deleted_sequence, tip_order_labels)
    stop(sprintf("deleted_sequence contains unknown tips: %s", paste(bad, collapse = ", ")))
  }

  max_delete_count <- length(deleted_sequence)
  if ((length(tip_order_labels) - max_delete_count) < min_retained_tip_count) {
    stop("fixed deleted_sequence violates min_retained_tip_count")
  }

  out <- vector("list", max_delete_count + 1L)
  for (step_id in 0:max_delete_count) {
    deleted_tip_labels <- if (step_id == 0L) character(0) else deleted_sequence[seq_len(step_id)]
    deleted_tip_label <- if (step_id == 0L) NA_character_ else deleted_sequence[step_id]
    retained_tip_labels <- setdiff(tip_order_labels, deleted_tip_labels)

    out[[step_id + 1L]] <- data.frame(
      step_id = step_id,
      schedule_type = schedule_name,
      start_rank = start_rank,
      deleted_tip_label = deleted_tip_label,
      retained_tip_count = length(retained_tip_labels),
      max_delete_count = max_delete_count,
      stringsAsFactors = FALSE
    )
    out[[step_id + 1L]]$deleted_tip_labels <- list(deleted_tip_labels)
    out[[step_id + 1L]]$retained_tip_labels <- list(retained_tip_labels)
  }

  schedule_tbl <- do.call(rbind, out)
  rownames(schedule_tbl) <- NULL
  schedule_tbl
}

flatten_schedule_table <- function(schedule_tbl) {
  out <- schedule_tbl
  if ("deleted_tip_labels" %in% colnames(out)) {
    out$deleted_tip_labels <- vapply(
      out$deleted_tip_labels,
      function(x) paste(x, collapse = ","),
      character(1)
    )
  }
  if ("retained_tip_labels" %in% colnames(out)) {
    out$retained_tip_labels <- vapply(
      out$retained_tip_labels,
      function(x) paste(x, collapse = ","),
      character(1)
    )
  }
  out
}

run_one_schedule <- function(full_tree2,
                             identity_tbl,
                             root_label,
                             tree_semantics,
                             tip_order_labels,
                             schedule_type,
                             fixed_deleted_sequence = NULL,
                             fixed_schedule_name = NULL,
                             start_rank,
                             delete_count,
                             delete_prop,
                             min_retained_tip_count,
                             max_delete_fraction,
                             run_id,
                             random_seed = NULL) {
  if (!is.null(fixed_deleted_sequence)) {
    schedule_tbl <- make_fixed_pruning_schedule(
      tip_order_labels = tip_order_labels,
      deleted_sequence = fixed_deleted_sequence,
      start_rank = start_rank,
      schedule_name = if (is.null(fixed_schedule_name)) "fixed_schedule" else fixed_schedule_name,
      min_retained_tip_count = min_retained_tip_count
    )
  } else {
    schedule_tbl <- make_pruning_schedule(
      tip_order_labels = tip_order_labels,
      schedule_type = schedule_type,
      start_rank = start_rank,
      delete_count = delete_count,
      delete_prop = delete_prop,
      min_retained_tip_count = min_retained_tip_count,
      max_delete_fraction = max_delete_fraction,
      random_seed = random_seed
    )
  }

  axis_ids <- as.character(identity_tbl$branch_id)
  state <- init_graph_state(identity_tbl, root_label, tree_semantics = tree_semantics)

  table_rows <- list()
  event_rows <- list()
  traj_rows <- list()
  cell_status_rows <- list()
  member_col_name <- if (identical(tree_semantics, "unrooted")) {
    "benchmark_unrooted_members"
  } else {
    "benchmark_rooted_members"
  }
  fusion_group_rows <- list()

  for (i in seq_len(nrow(schedule_tbl))) {
    step_id <- schedule_tbl$step_id[i]
    deleted_tip_label <- schedule_tbl$deleted_tip_label[i]

    if (!is.na(deleted_tip_label)) {
      res <- delete_tip_once(state, deleted_tip_label)
      state <- res$state
      event_rows[[length(event_rows) + 1L]] <- data.frame(
        run_id = run_id,
        step_id = step_id,
        deleted_tip = deleted_tip_label,
        deleted_tip_parent_label = res$event$deleted_tip_parent_label,
        local_e_term = res$event$local_e_term,
        local_e_sib = res$event$local_e_sib,
        local_e_up = res$event$local_e_up,
        local_sib_child_label = res$event$local_sib_child_label,
        local_sib_parent_label = res$event$local_sib_parent_label,
        local_up_child_label = res$event$local_up_child_label,
        local_up_parent_label = res$event$local_up_parent_label,
        stringsAsFactors = FALSE
      )
    }

    cls <- state_to_classification(state, identity_tbl)
    line <- matrix(cls$line[axis_ids], nrow = 1)
    colnames(line) <- axis_ids
    rownames(line) <- paste0(run_id, "_step", step_id)
    table_rows[[length(table_rows) + 1L]] <- line

    traj_core <- trajectory_summary_from_classification(cls, state)
    traj_rows[[length(traj_rows) + 1L]] <- data.frame(
      run_id = run_id,
      step_id = step_id,
      n_tips = schedule_tbl$retained_tip_count[i],
      deleted_tip = ifelse(is.na(deleted_tip_label), NA_character_, deleted_tip_label),
      traj_core,
      stringsAsFactors = FALSE
    )
    cell_status_rows[[length(cell_status_rows) + 1L]] <- classification_to_long_rows(
      classification = cls,
      run_id = run_id,
      step_id = step_id,
      semantics_label = tree_semantics
    )
    fg_rows <- state_to_fusion_group_rows(
      state = state,
      run_id = run_id,
      step_id = step_id,
      member_col_name = member_col_name
    )
    if (nrow(fg_rows) > 0L) {
      fusion_group_rows[[length(fusion_group_rows) + 1L]] <- fg_rows
    }

    covered <- unique(c(cls$na_struct, cls$na_fuse, names(cls$line)[!cls$line %in% c("NA_struct", "NA_fuse")]))
    if (!setequal(covered, axis_ids)) {
      stop(sprintf("Axis coverage failed at %s step %d", run_id, step_id))
    }
  }

  list(
    table0 = do.call(rbind, table_rows),
    oracle_events = if (length(event_rows) == 0L) data.frame() else do.call(rbind, event_rows),
    trajectory = do.call(rbind, traj_rows),
    oracle_cell_status_long = do.call(rbind, cell_status_rows),
    oracle_fusion_groups = if (length(fusion_group_rows) == 0L) {
      empty <- data.frame(
        gene_id = character(),
        step_id = integer(),
        merge_group_id = character(),
        expected_fused_length = character(),
        group_size = integer(),
        stringsAsFactors = FALSE
      )
      empty[[member_col_name]] <- character()
      empty
    } else {
      do.call(rbind, fusion_group_rows)
    },
    schedule_tbl = schedule_tbl,
    subtree_list = build_subtree_series(full_tree2, schedule_tbl)$subtree_list
  )
}

full_tree <- read_or_simulate_tree(TREE_FILE, N_TIPS_IF_SIM, SEED, TREE_SEMANTICS)
frozen <- freeze_full_tree_identity(full_tree, tree_semantics = TREE_SEMANTICS)
full_tree2 <- frozen$full_tree
identity_tbl <- frozen$identity_tbl
root_label <- frozen$root_label
tip_order_labels <- full_tree2$tip.label

scenario_preset <- scenario_presets[[SCENARIO_NAME]]
if (!is.null(scenario_preset)) {
  SCHEDULE_TYPE <- scenario_preset$schedule_name
  START_RANK <- scenario_preset$start_rank
  DELETE_COUNT <- length(scenario_preset$deleted_sequence)
  DELETE_PROP <- NULL
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (identical(SCHEDULE_TYPE, "random") && NREP_RANDOM > 0L) {
  all_tables <- list()
  all_events <- list()
  all_traj <- list()
  all_cell_status <- list()
  all_fusion_groups <- list()

  for (r in seq_len(NREP_RANDOM)) {
    run_id <- sprintf("rand%03d", r)
    rr <- run_one_schedule(
      full_tree2 = full_tree2,
      identity_tbl = identity_tbl,
      root_label = root_label,
      tree_semantics = TREE_SEMANTICS,
      tip_order_labels = tip_order_labels,
      schedule_type = "random",
      fixed_deleted_sequence = NULL,
      fixed_schedule_name = NULL,
      start_rank = 1L,
      delete_count = DELETE_COUNT,
      delete_prop = DELETE_PROP,
      min_retained_tip_count = MIN_RETAINED_TIPS,
      max_delete_fraction = MAX_DELETE_FRACTION,
      run_id = run_id,
      random_seed = SEED + r
    )
    all_tables[[r]] <- rr$table0
    all_events[[r]] <- rr$oracle_events
    all_traj[[r]] <- rr$trajectory
    all_cell_status[[r]] <- rr$oracle_cell_status_long
    all_fusion_groups[[r]] <- rr$oracle_fusion_groups
  }

  table_out <- do.call(rbind, all_tables)
  events_out <- do.call(rbind, all_events)
  traj_out <- do.call(rbind, all_traj)
  cell_status_out <- do.call(rbind, all_cell_status)
  fusion_groups_out <- do.call(rbind, all_fusion_groups)
} else {
  main_run <- run_one_schedule(
    full_tree2 = full_tree2,
    identity_tbl = identity_tbl,
    root_label = root_label,
    tree_semantics = TREE_SEMANTICS,
    tip_order_labels = tip_order_labels,
    schedule_type = SCHEDULE_TYPE,
    fixed_deleted_sequence = if (!is.null(scenario_preset)) scenario_preset$deleted_sequence else NULL,
    fixed_schedule_name = if (!is.null(scenario_preset)) scenario_preset$schedule_name else NULL,
    start_rank = START_RANK,
    delete_count = DELETE_COUNT,
    delete_prop = DELETE_PROP,
    min_retained_tip_count = MIN_RETAINED_TIPS,
    max_delete_fraction = MAX_DELETE_FRACTION,
    run_id = "main",
    random_seed = if (identical(SCHEDULE_TYPE, "random")) SEED else NULL
  )
  table_out <- main_run$table0
  events_out <- main_run$oracle_events
  traj_out <- main_run$trajectory
  cell_status_out <- main_run$oracle_cell_status_long
  fusion_groups_out <- main_run$oracle_fusion_groups
}

write_species_tree_file(
  full_tree2,
  paste0(OUT_PREFIX, ".species_tree.nwk")
)

if (exists("main_run")) {
  write_gene_tree_file(
    main_run$subtree_list,
    run_id = "main",
    outfile = paste0(OUT_PREFIX, ".gene_trees.nwk")
  )
}

write.table(
  table_out,
  file = paste0(OUT_PREFIX, ".table.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

write.table(
  table_out,
  file = file.path(OUT_DIR, "oracle_gene_by_original_branch_matrix.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

write.table(
  events_out,
  file = paste0(OUT_PREFIX, ".oracle_events.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  traj_out,
  file = paste0(OUT_PREFIX, ".trajectory.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  cell_status_out,
  file = file.path(OUT_DIR, "oracle_cell_status_long.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  fusion_groups_out,
  file = file.path(OUT_DIR, "oracle_fusion_groups.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

if (exists("main_run")) {
  write.table(
    flatten_schedule_table(main_run$schedule_tbl),
    file = paste0(OUT_PREFIX, ".schedule.tsv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

message("Done. Wrote: ",
        paste0(OUT_PREFIX, ".species_tree.nwk / ",
               OUT_PREFIX, ".gene_trees.nwk / ",
               OUT_PREFIX, ".table.txt / ",
               OUT_PREFIX, ".oracle_events.tsv / ",
               OUT_PREFIX, ".trajectory.tsv / ",
               file.path(OUT_DIR, "oracle_cell_status_long.tsv"), " / ",
               file.path(OUT_DIR, "oracle_fusion_groups.tsv"),
               if (exists("main_run")) paste0(" / ", OUT_PREFIX, ".schedule.tsv") else ""))

plot_script <- file.path(script_dir, "plot_fulltree_collapse.R")
if (!file.exists(plot_script)) {
  stop(sprintf("Cannot find plotting script: %s", plot_script))
}
source(plot_script, local = FALSE)
