#ifndef SPLITALIGNER_FAST_C_API_H
#define SPLITALIGNER_FAST_C_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#define SAF_ABI_VERSION 2u
#define SAF_NONE ((size_t)-1)
typedef struct saf_reference saf_reference;
typedef struct saf_coverage saf_coverage;
typedef struct saf_mapping saf_mapping;
typedef struct saf_support saf_support;
typedef enum saf_status { SAF_OK=0, SAF_INVALID_ARGUMENT=1, SAF_CAPACITY=2, SAF_NO_MEMORY=3, SAF_NUMERIC_OVERFLOW=4, SAF_COUNT_OVERFLOW=5, SAF_INTERNAL_ERROR=6 } saf_status;
typedef struct saf_reference_info { size_t n_taxa, n_words, n_axes, n_nodes; } saf_reference_info;
typedef struct saf_coverage_info { size_t n_taxa, n_words, n_axes, n_projections, n_observed; } saf_coverage_info;
typedef struct saf_mapping_info { size_t n_nodes, n_axes, n_projections, n_gene_splits, n_reference_absent, n_gene_unmatched; } saf_mapping_info;
typedef struct saf_split_info { size_t side_count, other_count; uint8_t terminal; } saf_split_info;
typedef struct saf_value_buffer { double *values; uint8_t *valid; size_t capacity; } saf_value_buffer;
typedef struct saf_support_count { uint64_t n_valid_topologies, n_struct, n_fused, n_dec, n_rec, n_not_recovered; } saf_support_count;
const char *saf_version(void);
const char *saf_last_error(void);
/* parent[node]=-1 for root; taxon[node]=-1 for internal; inputs need not be ordered. */
saf_status saf_reference_create(size_t n_taxa, const char *universe_id, const int64_t *parent, const int64_t *taxon, size_t n_nodes, saf_reference **out);
void saf_reference_destroy(saf_reference *h);
saf_status saf_reference_get_info(const saf_reference *h, saf_reference_info *out);
saf_status saf_reference_axis(const saf_reference *h, size_t axis, uint64_t *words, size_t capacity, size_t *out_count, saf_split_info *info);
saf_status saf_reference_sources(const saf_reference *h, size_t axis, size_t *nodes, size_t capacity, size_t *out_count);
saf_status saf_reference_universe(const saf_reference *h, char *bytes, size_t capacity, size_t *out_count);
saf_status saf_coverage_create(const saf_reference *reference, const uint64_t *mask, size_t n_words, saf_coverage **out);
void saf_coverage_destroy(saf_coverage *h);
saf_status saf_coverage_get_info(const saf_coverage *h, saf_coverage_info *out);
saf_status saf_coverage_mask(const saf_coverage *h, uint64_t *words, size_t capacity, size_t *out_count);
saf_status saf_coverage_projection(const saf_coverage *h, size_t projection, uint64_t *words, size_t capacity, size_t *out_count, saf_split_info *info);
saf_status saf_coverage_members(const saf_coverage *h, size_t projection, size_t *axes, size_t capacity, size_t *out_count);
/* Fixed vectors: all three arrays mandatory when n_axes>0; validate capacity before copying. */
saf_status saf_coverage_axes(const saf_coverage *h, size_t *axis_to_projection, uint8_t *structural, uint8_t *eligible, size_t capacity);
saf_status saf_mapping_create(const saf_coverage *coverage, const int64_t *parent, const int64_t *taxon, size_t n_nodes, saf_mapping **out);
void saf_mapping_destroy(saf_mapping *h);
saf_status saf_mapping_get_info(const saf_mapping *h, saf_mapping_info *out);
saf_status saf_mapping_gene_split(const saf_mapping *h, size_t gene_edge, uint64_t *words, size_t capacity, size_t *out_count, saf_split_info *info);
saf_status saf_mapping_gene_sources(const saf_mapping *h, size_t gene_edge, size_t *nodes, size_t capacity, size_t *out_count);
saf_status saf_mapping_projection_to_gene(const saf_mapping *h, size_t *indices, size_t capacity, size_t *out_count);
saf_status saf_mapping_gene_to_projection(const saf_mapping *h, size_t *indices, size_t capacity, size_t *out_count);
saf_status saf_mapping_reference_absent(const saf_mapping *h, size_t *indices, size_t capacity, size_t *out_count);
saf_status saf_mapping_gene_unmatched(const saf_mapping *h, size_t *indices, size_t capacity, size_t *out_count);
/* raw_states has n_axes bytes S/F/D/U; recovered has n_projections bytes 0/1. */
saf_status saf_mapping_states(const saf_mapping *h, char *raw_states, size_t capacity);
saf_status saf_mapping_recovered(const saf_mapping *h, uint8_t *recovered, size_t capacity);
/* Skip a family with {NULL,NULL,0}; requested families require both pointers. */
saf_status saf_mapping_apply(const saf_mapping *h, const double *node_lengths, const uint8_t *length_valid, size_t n_nodes, saf_value_buffer primitive, saf_value_buffer projection, saf_value_buffer gene_edge);
saf_status saf_support_create(const saf_reference *reference, const char *unit_id, const char *cohort_id, saf_support **out);
void saf_support_destroy(saf_support *h);
saf_status saf_support_update(saf_support *h, const saf_mapping *observation);
saf_status saf_support_merge(saf_support *h, const saf_support *disjoint_shard);
saf_status saf_support_read(const saf_support *h, saf_support_count *counts, size_t capacity, size_t *out_count);
#ifdef __cplusplus
}
#endif
#endif
