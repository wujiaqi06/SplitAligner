#include "splitaligner/fast_c_api.h"
#include <math.h>
#include <stdio.h>
#define REQUIRE(x) do { if (!(x)) { fprintf(stderr,"C ABI failure line %d: %s\n",__LINE__,saf_last_error()); return 1; } } while(0)
int main(void) {
 int64_t parent[3]={-1,0,0},taxon[3]={-1,0,1};
 saf_reference *r=NULL;saf_coverage *c=NULL;saf_mapping *m=NULL;saf_support *s=NULL;
 uint64_t mask=3;double lengths[3]={NAN,0.25,0.5},value[2]={-9,-77};uint8_t valid[2]={0,99};
 saf_value_buffer request={value,valid,2},skip={NULL,NULL,0};size_t n=0,source[2];saf_support_count count;
 REQUIRE(saf_reference_create(2,"C-smoke:0=a;1=b",parent,taxon,3,&r)==SAF_OK);
 REQUIRE(saf_coverage_create(r,&mask,1,&c)==SAF_OK);
 REQUIRE(saf_mapping_create(c,parent,taxon,3,&m)==SAF_OK);
 REQUIRE(saf_support_create(r,"gene","c-smoke",&s)==SAF_OK);
 saf_reference_destroy(r);saf_coverage_destroy(c);
 REQUIRE(saf_mapping_gene_sources(m,0,NULL,0,&n)==SAF_OK && n==2);
 REQUIRE(saf_mapping_gene_sources(m,0,source,2,&n)==SAF_OK && source[0]==1 && source[1]==2);
 REQUIRE(saf_mapping_apply(m,lengths,NULL,3,request,skip,skip)==SAF_OK);
 REQUIRE(value[0]==0.75 && valid[0]==1 && value[1]==-77 && valid[1]==99);
 REQUIRE(saf_support_update(s,m)==SAF_OK);
 REQUIRE(saf_support_read(s,&count,1,&n)==SAF_OK && n==1 && count.n_dec==1 && count.n_rec==1);
 saf_mapping_destroy(m);saf_support_destroy(s);saf_reference_destroy(NULL);saf_coverage_destroy(NULL);saf_mapping_destroy(NULL);saf_support_destroy(NULL);
 puts("native_c_smoke: PASS");return 0;
}
