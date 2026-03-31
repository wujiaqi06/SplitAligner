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

# -----------------------------
# User parameters (edit here)
# -----------------------------

SEED <- 42
N_TIPS_IF_SIM <- 10

# set TREE_FILE to a Newick path; otherwise a random tree is simulated
TREE_FILE <- file.path(base_dir, "inputs", "example_tree.nwk")

# schedule types:
#   global_ordered / local_forward / local_backward / random
SCHEDULE_TYPE <- "global_ordered"
START_RANK <- 1L
DELETE_COUNT <- NULL
DELETE_PROP <- NULL
MIN_RETAINED_TIPS <- 3L
# Benchmark V1 packaged default.
# This packaged demo uses 0.7 to make collapse / fusion progression more visible.
MAX_DELETE_FRACTION <- 0.7

# if > 0 and SCHEDULE_TYPE == "random", multiple replicates are run
NREP_RANDOM <- 0L

OUT_PREFIX <- file.path(base_dir, "outputs", "benchmark")

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

prepare_benchmark_tree <- function(tr) {
  if (is.null(tr$edge.length)) {
    stop("Benchmark tree must have branch lengths")
  }
  if (anyNA(tr$edge.length)) {
    stop("Benchmark tree contains missing branch lengths")
  }
  if (is.rooted(tr)) {
    tr <- unroot(tr)
  }
  tr
}

read_or_simulate_tree <- function(tree_file, n_tips, seed) {
  set.seed(seed)
  if (!is.null(tree_file) && file.exists(tree_file)) {
    tr <- read.tree(tree_file)
    if (inherits(tr, "multiPhylo")) stop("TREE_FILE must contain exactly one tree")
    return(prepare_benchmark_tree(tr))
  }
  prepare_benchmark_tree(rtree(n_tips))
}

run_one_schedule <- function(full_tree2,
                             identity_tbl,
                             root_label,
                             tip_order_labels,
                             schedule_type,
                             start_rank,
                             delete_count,
                             delete_prop,
                             min_retained_tip_count,
                             max_delete_fraction,
                             run_id,
                             random_seed = NULL) {
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

  axis_ids <- as.character(identity_tbl$branch_id)
  state <- init_graph_state(identity_tbl, root_label)

  table_rows <- list()
  event_rows <- list()
  traj_rows <- list()

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

    covered <- unique(c(cls$na_struct, cls$na_fuse, names(cls$line)[!cls$line %in% c("NA_struct", "NA_fuse")]))
    if (!setequal(covered, axis_ids)) {
      stop(sprintf("Axis coverage failed at %s step %d", run_id, step_id))
    }
  }

  list(
    table0 = do.call(rbind, table_rows),
    oracle_events = if (length(event_rows) == 0L) data.frame() else do.call(rbind, event_rows),
    trajectory = do.call(rbind, traj_rows),
    schedule_tbl = schedule_tbl,
    subtree_list = build_subtree_series(full_tree2, schedule_tbl)$subtree_list
  )
}

full_tree <- read_or_simulate_tree(TREE_FILE, N_TIPS_IF_SIM, SEED)
frozen <- freeze_full_tree_identity(full_tree)
full_tree2 <- frozen$full_tree
identity_tbl <- frozen$identity_tbl
root_label <- frozen$root_label
tip_order_labels <- full_tree2$tip.label

if (identical(SCHEDULE_TYPE, "random") && NREP_RANDOM > 0L) {
  all_tables <- list()
  all_events <- list()
  all_traj <- list()

  for (r in seq_len(NREP_RANDOM)) {
    run_id <- sprintf("rand%03d", r)
    rr <- run_one_schedule(
      full_tree2 = full_tree2,
      identity_tbl = identity_tbl,
      root_label = root_label,
      tip_order_labels = tip_order_labels,
      schedule_type = "random",
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
  }

  table_out <- do.call(rbind, all_tables)
  events_out <- do.call(rbind, all_events)
  traj_out <- do.call(rbind, all_traj)
} else {
  main_run <- run_one_schedule(
    full_tree2 = full_tree2,
    identity_tbl = identity_tbl,
    root_label = root_label,
    tip_order_labels = tip_order_labels,
    schedule_type = SCHEDULE_TYPE,
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

message("Done. Wrote: ",
        paste0(OUT_PREFIX, ".species_tree.nwk / ",
               OUT_PREFIX, ".gene_trees.nwk / ",
               OUT_PREFIX, ".table.txt / ",
               OUT_PREFIX, ".oracle_events.tsv / ",
               OUT_PREFIX, ".trajectory.tsv"))

plot_script <- file.path(script_dir, "plot_fulltree_collapse.R")
if (!file.exists(plot_script)) {
  stop(sprintf("Cannot find plotting script: %s", plot_script))
}
source(plot_script, local = FALSE)
