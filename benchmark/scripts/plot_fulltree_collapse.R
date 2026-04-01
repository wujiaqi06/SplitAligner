suppressPackageStartupMessages(library(ape))

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
tree_semantics <- if (exists("TREE_SEMANTICS", inherits = TRUE)) {
  get("TREE_SEMANTICS", inherits = TRUE)
} else {
  "unrooted"
}
out_dir <- file.path(base_dir, "outputs", tree_semantics)
out_prefix <- if (exists("OUT_PREFIX", inherits = TRUE)) {
  get("OUT_PREFIX", inherits = TRUE)
} else {
  file.path(out_dir, "benchmark")
}

TREE_FILE <- paste0(out_prefix, ".species_tree.nwk")
TABLE_TXT <- paste0(out_prefix, ".table.txt")
EVENTS_TSV <- paste0(out_prefix, ".oracle_events.tsv")

OUT_PDF <- file.path(out_dir, "benchmark.fulltree_collapse.pdf")

tr0 <- read.tree(TREE_FILE)

tab <- read.table(TABLE_TXT, header=TRUE, row.names=1, sep="\t",
                  check.names=FALSE, stringsAsFactors=FALSE)

events <- read.delim(EVENTS_TSV, stringsAsFactors=FALSE, check.names=FALSE)
events$key <- paste0(events$run_id, "_step", events$step_id)  # e.g. main_step1

# --- helpers ---------------------------------------------------------------

node_num_from_label <- function(tr, lab) {
  if (is.na(lab) || lab == "" || lab == "NA") return(NA_integer_)
  ti <- match(lab, tr$tip.label)
  if (!is.na(ti)) return(as.integer(ti))
  if (!is.null(tr$node.label)) {
    ni <- match(lab, tr$node.label)
    if (!is.na(ni)) return(as.integer(length(tr$tip.label) + ni))
  }
  NA_integer_
}

edge_index_from_child_label <- function(tr, child_lab) {
  ch <- node_num_from_label(tr, child_lab)
  if (is.na(ch)) return(NA_integer_)
  idx <- which(tr$edge[,2] == ch)
  if (length(idx) == 0) return(NA_integer_)
  as.integer(idx[1])
}

get_state <- function(x) {
  if (is.na(x)) return("unknown")
  if (x == "NA_struct") return("na_struct")
  if (x == "NA_fuse") return("na_fuse")
  suppressWarnings({
    v <- as.numeric(x)
    if (!is.na(v)) return("present")
  })
  "unknown"
}

# --- build mapping edge_id -> edge index on T0 ----------------------------

edge_ids <- colnames(tab)
edge_idx <- sapply(edge_ids, function(e) edge_index_from_child_label(tr0, e))
valid <- which(!is.na(edge_idx))
edge_ids_v <- edge_ids[valid]
edge_idx_v <- edge_idx[valid]

# --- plotting: cumulative collapse on full tree ---------------------------

pdf(OUT_PDF, width=11, height=8.5)

cum_struct <- character(0)
cum_fuse   <- character(0)
deleted_so_far <- character(0)

for (r in rownames(tab)) {
  
  # table-driven cumulative update
  row <- tab[r, edge_ids_v, drop=TRUE]
  states <- sapply(row, get_state)
  
  newly_struct <- edge_ids_v[states == "na_struct"]
  newly_fuse   <- edge_ids_v[states == "na_fuse"]
  
  cum_struct <- union(cum_struct, newly_struct)
  cum_fuse   <- union(cum_fuse, newly_fuse)
  cum_fuse <- setdiff(cum_fuse, cum_struct)
  
  # event row for this step
  ev_i <- events[events$key == r, ]
  
  # update deleted tips shown in title
  if (nrow(ev_i) == 1) {
    deleted_so_far <- unique(c(deleted_so_far, as.character(ev_i$deleted_tip)))
  }
  
  # root symmetry patch
  if (identical(tree_semantics, "rooted") && nrow(ev_i) == 1) {
    if (!is.na(ev_i$deleted_tip_parent_label) && ev_i$deleted_tip_parent_label == "__ROOT__") {
      extra <- as.character(ev_i$local_e_sib)
      if (!is.na(extra) && extra != "" && extra != "NA") {
        cum_struct <- union(cum_struct, extra)
        cum_fuse <- setdiff(cum_fuse, cum_struct)
      }
    }
  }
  
  # edge styles
  nedge <- nrow(tr0$edge)
  e_col <- rep("black", nedge)
  e_lwd <- rep(1.2, nedge)
  e_lty <- rep(1, nedge)
  
  if (length(cum_fuse) > 0) {
    ii <- edge_idx_v[match(cum_fuse, edge_ids_v)]
    ii <- ii[!is.na(ii)]
    e_col[ii] <- "blue"
    e_lwd[ii] <- 3.0
  }
  if (length(cum_struct) > 0) {
    ii <- edge_idx_v[match(cum_struct, edge_ids_v)]
    ii <- ii[!is.na(ii)]
    e_col[ii] <- "red"
    e_lwd[ii] <- 3.2
  }
  
  plot(tr0,
       show.tip.label = TRUE,
       edge.color = e_col,
       edge.width = e_lwd,
       edge.lty = e_lty,
       cex = 0.9)
  
  title(main = sprintf("Collapse on full species tree: %s | deleted={%s}",
                       r, paste(deleted_so_far, collapse = ",")))
  
  legend("topleft", bty="n", cex=0.9,
         legend = c(sprintf("NA_struct (cumulative): %d", length(cum_struct)),
                    sprintf("NA_fuse (cumulative): %d", length(cum_fuse)),
                    "unaffected"),
         text.col = c("red","blue","black"),
         lwd = c(3.2,3.0,1.2))
}

dev.off()

cat("Done. PDF written to:", OUT_PDF, "\n")
