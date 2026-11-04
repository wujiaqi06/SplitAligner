#include "splitaligner/fast_core.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <new>
#include <utility>
struct saf_reference {std::shared_ptr<const saf::Reference> p;};
struct saf_coverage {std::shared_ptr<const saf::Coverage> p;};
struct saf_mapping {std::shared_ptr<const saf::Mapping> p;};
struct saf_support {std::unique_ptr<saf::Support> p;};
namespace {
thread_local char error_text[512]={0};
void error(const char*m)noexcept{size_t i=0;if(m)for(;i+1<sizeof(error_text)&&m[i];++i)error_text[i]=m[i];error_text[i]=0;}
[[noreturn]] void invalid(const char*m){throw saf::Error(SAF_INVALID_ARGUMENT,m);}
template<class F>saf_status guard(F f)noexcept{try{f();error_text[0]=0;return SAF_OK;}catch(const saf::Error&e){error(e.what());return e.code;}catch(const std::bad_alloc&){error("allocation failed");return SAF_NO_MEMORY;}catch(const std::length_error&e){error(e.what());return SAF_INVALID_ARGUMENT;}catch(const std::exception&e){error(e.what());return SAF_INTERNAL_ERROR;}catch(...){error("unknown C++ exception");return SAF_INTERNAL_ERROR;}}
template<class T>const auto& need(const T*h){if(!h||!h->p)invalid("null handle");return *h->p;}
template<class T>void query(const T*source,size_t n,T*out,size_t cap,size_t*count){if(!count)invalid("out_count is required");*count=n;if(!out){if(cap)invalid("NULL query buffer needs capacity zero");return;}if(cap<n)throw saf::Error(SAF_CAPACITY,"query capacity too small");if(n)std::copy(source,source+n,out);}
template<class T>void query(const std::vector<T>&v,T*out,size_t cap,size_t*count){query(v.data(),v.size(),out,cap,count);}
void fixed(const void*out,size_t cap,size_t n){if(n&&!out)invalid("fixed output buffer required");if(cap<n)throw saf::Error(SAF_CAPACITY,"fixed output capacity too small");if(!out&&cap)invalid("NULL fixed output requires zero capacity");}
const saf::Split& axis(const saf::Reference&r,size_t i){if(i>=r.axes().size())invalid("reference axis out of range");return r.axes()[i];}
const saf::Projection& projection(const saf::Coverage&c,size_t i){if(i>=c.projections().size())invalid("projection out of range");return c.projections()[i];}
const saf::Split& edge(const saf::Mapping&m,size_t i){if(i>=m.gene_splits().size())invalid("gene edge out of range");return m.gene_splits()[i];}
void split_query(const saf::Split&s,uint64_t*out,size_t cap,size_t*count,saf_split_info*info){query(s.side,out,cap,count);if(info)*info={s.side_count,s.other_count,static_cast<uint8_t>(s.terminal)};}
}
extern "C" {
const char*saf_version(void){return "SplitAlignerFast FAST002 ABI2 topology-support-v1 finite-double-v1";}
const char*saf_last_error(void){return error_text;}
saf_status saf_reference_create(size_t nt,const char*universe,const int64_t*p,const int64_t*t,size_t nn,saf_reference**out){if(out)*out=nullptr;return guard([&]{if(!out)invalid("out_handle required");auto plan=saf::Reference::create(nt,universe,p,t,nn);*out=new saf_reference{std::move(plan)};});}
void saf_reference_destroy(saf_reference*h){guard([&]{delete h;});}
saf_status saf_reference_get_info(const saf_reference*h,saf_reference_info*out){return guard([&]{const auto&r=need(h);if(!out)invalid("info output required");*out={r.n_taxa(),r.n_words(),r.axes().size(),r.n_nodes()};});}
saf_status saf_reference_axis(const saf_reference*h,size_t i,uint64_t*out,size_t cap,size_t*count,saf_split_info*info){return guard([&]{split_query(axis(need(h),i),out,cap,count,info);});}
saf_status saf_reference_sources(const saf_reference*h,size_t i,size_t*out,size_t cap,size_t*count){return guard([&]{query(axis(need(h),i).sources,out,cap,count);});}
saf_status saf_reference_universe(const saf_reference*h,char*out,size_t cap,size_t*count){return guard([&]{const auto&s=need(h).universe();query(s.c_str(),s.size()+1,out,cap,count);});}
saf_status saf_coverage_create(const saf_reference*r,const uint64_t*mask,size_t nw,saf_coverage**out){if(out)*out=nullptr;return guard([&]{if(!out)invalid("out_handle required");need(r);auto plan=saf::Coverage::create(r->p,mask,nw);*out=new saf_coverage{std::move(plan)};});}
void saf_coverage_destroy(saf_coverage*h){guard([&]{delete h;});}
saf_status saf_coverage_get_info(const saf_coverage*h,saf_coverage_info*out){return guard([&]{const auto&c=need(h);if(!out)invalid("info output required");*out={c.reference()->n_taxa(),c.mask().size(),c.reference()->axes().size(),c.projections().size(),c.observed()};});}
saf_status saf_coverage_mask(const saf_coverage*h,uint64_t*out,size_t cap,size_t*count){return guard([&]{query(need(h).mask(),out,cap,count);});}
saf_status saf_coverage_projection(const saf_coverage*h,size_t i,uint64_t*out,size_t cap,size_t*count,saf_split_info*info){return guard([&]{split_query(projection(need(h),i).split,out,cap,count,info);});}
saf_status saf_coverage_members(const saf_coverage*h,size_t i,size_t*out,size_t cap,size_t*count){return guard([&]{query(projection(need(h),i).members,out,cap,count);});}
saf_status saf_coverage_axes(const saf_coverage*h,size_t*indices,uint8_t*structural,uint8_t*eligible,size_t cap){return guard([&]{const auto&c=need(h);size_t n=c.axis_to_projection().size();fixed(indices,cap,n);fixed(structural,cap,n);fixed(eligible,cap,n);std::copy(c.axis_to_projection().begin(),c.axis_to_projection().end(),indices);std::copy(c.structural().begin(),c.structural().end(),structural);std::copy(c.eligible().begin(),c.eligible().end(),eligible);});}
saf_status saf_mapping_create(const saf_coverage*c,const int64_t*p,const int64_t*t,size_t nn,saf_mapping**out){if(out)*out=nullptr;return guard([&]{if(!out)invalid("out_handle required");need(c);auto plan=saf::Mapping::create(c->p,p,t,nn);*out=new saf_mapping{std::move(plan)};});}
void saf_mapping_destroy(saf_mapping*h){guard([&]{delete h;});}
saf_status saf_mapping_get_info(const saf_mapping*h,saf_mapping_info*out){return guard([&]{const auto&m=need(h);if(!out)invalid("info output required");*out={m.n_nodes(),m.states().size(),m.projection_to_gene().size(),m.gene_splits().size(),m.reference_absent().size(),m.gene_unmatched().size()};});}
saf_status saf_mapping_gene_split(const saf_mapping*h,size_t i,uint64_t*out,size_t cap,size_t*count,saf_split_info*info){return guard([&]{split_query(edge(need(h),i),out,cap,count,info);});}
saf_status saf_mapping_gene_sources(const saf_mapping*h,size_t i,size_t*out,size_t cap,size_t*count){return guard([&]{query(edge(need(h),i).sources,out,cap,count);});}
saf_status saf_mapping_projection_to_gene(const saf_mapping*h,size_t*out,size_t cap,size_t*count){return guard([&]{query(need(h).projection_to_gene(),out,cap,count);});}
saf_status saf_mapping_gene_to_projection(const saf_mapping*h,size_t*out,size_t cap,size_t*count){return guard([&]{query(need(h).gene_to_projection(),out,cap,count);});}
saf_status saf_mapping_reference_absent(const saf_mapping*h,size_t*out,size_t cap,size_t*count){return guard([&]{query(need(h).reference_absent(),out,cap,count);});}
saf_status saf_mapping_gene_unmatched(const saf_mapping*h,size_t*out,size_t cap,size_t*count){return guard([&]{query(need(h).gene_unmatched(),out,cap,count);});}
saf_status saf_mapping_states(const saf_mapping*h,char*out,size_t cap){return guard([&]{const auto&v=need(h).states();fixed(out,cap,v.size());if(!v.empty())std::copy(v.begin(),v.end(),out);});}
saf_status saf_mapping_recovered(const saf_mapping*h,uint8_t*out,size_t cap){return guard([&]{const auto&v=need(h).recovered();fixed(out,cap,v.size());if(!v.empty())std::copy(v.begin(),v.end(),out);});}
saf_status saf_mapping_apply(const saf_mapping*h,const double*values,const uint8_t*valid,size_t nn,saf_value_buffer primitive,saf_value_buffer projection,saf_value_buffer gene){return guard([&]{need(h).apply(values,valid,nn,primitive,projection,gene);});}
saf_status saf_support_create(const saf_reference*r,const char*unit,const char*cohort,saf_support**out){if(out)*out=nullptr;return guard([&]{if(!out)invalid("out_handle required");need(r);auto plan=saf::Support::create(r->p,unit,cohort);*out=new saf_support{std::move(plan)};});}
void saf_support_destroy(saf_support*h){guard([&]{delete h;});}
saf_status saf_support_update(saf_support*h,const saf_mapping*m){return guard([&]{need(h);h->p->update(need(m));});}
saf_status saf_support_merge(saf_support*h,const saf_support*s){return guard([&]{need(h);h->p->merge(need(s));});}
saf_status saf_support_read(const saf_support*h,saf_support_count*out,size_t cap,size_t*count){return guard([&]{query(need(h).counts(),out,cap,count);});}
}
