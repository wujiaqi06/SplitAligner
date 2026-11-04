#include "splitaligner/fast_c_api.h"
#include "splitaligner/fast_core.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

static_assert(!std::is_copy_constructible<saf::Reference>::value, "Reference must not be copied");
static_assert(!std::is_move_constructible<saf::Reference>::value, "Reference must not be moved");
static_assert(!std::is_copy_assignable<saf::Reference>::value, "Reference must not be copy-assigned");
static_assert(!std::is_move_assignable<saf::Reference>::value, "Reference must not be move-assigned");
static_assert(!std::is_copy_constructible<saf::Coverage>::value, "Coverage must not be copied");
static_assert(!std::is_move_constructible<saf::Coverage>::value, "Coverage must not be moved");
static_assert(!std::is_copy_assignable<saf::Coverage>::value, "Coverage must not be copy-assigned");
static_assert(!std::is_move_assignable<saf::Coverage>::value, "Coverage must not be move-assigned");
static_assert(!std::is_copy_constructible<saf::Mapping>::value, "Mapping must not be copied");
static_assert(!std::is_move_constructible<saf::Mapping>::value, "Mapping must not be moved");
static_assert(!std::is_copy_assignable<saf::Mapping>::value, "Mapping must not be copy-assigned");
static_assert(!std::is_move_assignable<saf::Mapping>::value, "Mapping must not be move-assigned");

#ifdef SAF_COMPILE_FAIL_WRITABLE_COPY
// Negative compilation probe: this formerly ordinary C++ copy could retain a
// writable alias to the Reference held by a Coverage and later replace its state.
// The repaired header must reject both make_shared's copy and the assignment.
void prior_writable_copy_probe(std::shared_ptr<const saf::Reference> original,
                               std::shared_ptr<const saf::Reference> replacement) {
 auto writable = std::make_shared<saf::Reference>(*original);
 uint64_t mask = 3;
 auto coverage = saf::Coverage::create(writable, &mask, 1);
 *writable = *replacement;
 (void)coverage;
}
#endif

static std::atomic<size_t> allocation_count{0};
static std::atomic<bool> track_allocations{false};
static std::atomic<bool> fail_next_allocation{false};
void* operator new(size_t n){if(fail_next_allocation.exchange(false))throw std::bad_alloc();if(track_allocations.load())++allocation_count;void*p=std::malloc(n?n:1);if(!p)throw std::bad_alloc();return p;}
void* operator new[](size_t n){return ::operator new(n);}
void operator delete(void*p)noexcept{std::free(p);}
void operator delete[](void*p)noexcept{std::free(p);}
void operator delete(void*p,size_t)noexcept{std::free(p);}
void operator delete[](void*p,size_t)noexcept{std::free(p);}
#define CHECK(x) do {if(!(x)){std::fprintf(stderr,"native failure line %d: %s; %s\n",__LINE__,#x,saf_last_error());std::abort();}}while(0)
#define OK(x) CHECK((x)==SAF_OK)
const saf_value_buffer skip{nullptr,nullptr,0};
struct Fixture {
 saf_reference*r=nullptr;saf_coverage*c=nullptr;saf_mapping*m=nullptr;
 ~Fixture(){saf_mapping_destroy(m);saf_coverage_destroy(c);saf_reference_destroy(r);}
};
void two(Fixture&f){int64_t p[]={-1,0,0},t[]={-1,0,1};uint64_t mask=3;OK(saf_reference_create(2,"two:0=a;1=b",p,t,3,&f.r));OK(saf_coverage_create(f.r,&mask,1,&f.c));OK(saf_mapping_create(f.c,p,t,3,&f.m));}
void numeric_and_queries(){
 Fixture f;two(f);size_t count=0;uint64_t word=999;saf_split_info si{};
 OK(saf_reference_axis(f.r,0,nullptr,0,&count,&si));CHECK(count==1&&si.terminal&&si.side_count==1&&si.other_count==1);
 CHECK(saf_reference_axis(f.r,0,&word,0,&count,&si)==SAF_CAPACITY&&word==999&&count==1);
 CHECK(saf_reference_axis(f.r,0,&word,1,nullptr,&si)==SAF_INVALID_ARGUMENT);
 CHECK(saf_reference_axis(f.r,0,nullptr,1,&count,&si)==SAF_INVALID_ARGUMENT);
 size_t sources[3]={99,99,99};CHECK(saf_mapping_gene_sources(f.m,0,sources,1,&count)==SAF_CAPACITY&&count==2&&sources[0]==99);
 OK(saf_mapping_gene_sources(f.m,0,sources,3,&count));CHECK(sources[0]==1&&sources[1]==2&&sources[2]==99);
 double input[3]={NAN,0.25,0.5},o[2]={99,77};uint8_t av[2]={99,77};saf_value_buffer out{o,av,2};
 OK(saf_mapping_apply(f.m,nullptr,nullptr,0,skip,skip,skip));
 CHECK(saf_mapping_apply(f.m,nullptr,nullptr,0,{nullptr,nullptr,1},skip,skip)==SAF_INVALID_ARGUMENT);
 CHECK(saf_mapping_apply(f.m,input,nullptr,3,{o,av,0},skip,skip)==SAF_CAPACITY&&o[0]==99);
 CHECK(saf_mapping_apply(f.m,input,nullptr,2,out,skip,skip)==SAF_INVALID_ARGUMENT&&o[0]==99);
 CHECK(saf_mapping_apply(f.m,input,nullptr,3,{input,av,2},skip,skip)==SAF_INVALID_ARGUMENT);
 CHECK(saf_mapping_apply(f.m,input,nullptr,3,out,out,skip)==SAF_INVALID_ARGUMENT);
 OK(saf_mapping_apply(f.m,input,nullptr,3,out,skip,skip));CHECK(o[0]==0.75&&av[0]&&o[1]==77&&av[1]==77);
 double qo=0,go=0;uint8_t qv=0,gv=0;
 allocation_count=0;track_allocations=true;
 for(int i=0;i<10000;++i){if(saf_mapping_apply(f.m,input,nullptr,3,out,{&qo,&qv,1},{&go,&gv,1})!=SAF_OK)std::abort();}
 track_allocations=false;CHECK(allocation_count==0);
 input[1]=-0.25;input[2]=0;OK(saf_mapping_apply(f.m,input,nullptr,3,skip,skip,out));CHECK(o[0]==-0.25&&av[0]);
 input[1]=std::numeric_limits<double>::denorm_min();OK(saf_mapping_apply(f.m,input,nullptr,3,skip,out,skip));CHECK(o[0]==input[1]&&av[0]);
 input[1]=1e308;input[2]=1e308;CHECK(saf_mapping_apply(f.m,input,nullptr,3,out,skip,skip)==SAF_NUMERIC_OVERFLOW);
 uint8_t available[3]={0,1,0};OK(saf_mapping_apply(f.m,input,available,3,out,skip,skip));CHECK(!av[0]&&std::isnan(o[0]));
 available[2]=1;input[1]=NAN;OK(saf_mapping_apply(f.m,input,available,3,out,skip,skip));CHECK(!av[0]&&std::isnan(o[0]));
 input[1]=INFINITY;OK(saf_mapping_apply(f.m,input,nullptr,3,out,skip,skip));CHECK(!av[0]);
 input[1]=0.25;input[2]=0.5;OK(saf_mapping_apply(f.m,input,nullptr,3,out,skip,skip));CHECK(o[0]==0.75);
 int64_t reordered_p[]={2,2,-1},reordered_t[]={0,1,-1};saf_mapping*reordered=nullptr;
 OK(saf_mapping_create(f.c,reordered_p,reordered_t,3,&reordered));double reordered_lengths[]={0.25,0.5,NAN};
 OK(saf_mapping_apply(reordered,reordered_lengths,nullptr,3,out,skip,skip));CHECK(o[0]==0.75);saf_mapping_destroy(reordered);
 saf_reference_destroy(f.r);f.r=nullptr;saf_coverage_destroy(f.c);f.c=nullptr;
 OK(saf_mapping_apply(f.m,input,nullptr,3,out,skip,skip));CHECK(o[0]==0.75);
 std::thread a([&]{for(int i=0;i<1000;++i){double local;uint8_t v;OK(saf_mapping_apply(f.m,input,nullptr,3,{&local,&v,1},skip,skip));CHECK(local==0.75);}});
 std::thread b([&]{for(int i=0;i<1000;++i){double local;uint8_t v;OK(saf_mapping_apply(f.m,input,nullptr,3,{&local,&v,1},skip,skip));CHECK(local==0.75);}});a.join();b.join();
}
void fusion_and_topology(){
 Fixture f;int64_t rp[]={-1,0,1,1,0,0,5,5},rt[]={-1,-1,0,1,4,-1,2,3};uint64_t mask=15;
 OK(saf_reference_create(5,"five:A,B,C,D,E",rp,rt,8,&f.r));OK(saf_coverage_create(f.r,&mask,1,&f.c));
 int64_t gp[]={-1,0,1,1,0,4,4},gt[]={-1,-1,0,1,-1,2,3};OK(saf_mapping_create(f.c,gp,gt,7,&f.m));
 saf_mapping_info mi;OK(saf_mapping_get_info(f.m,&mi));CHECK(mi.n_axes==7&&mi.n_projections==5&&mi.n_gene_splits==5&&mi.n_reference_absent==0&&mi.n_gene_unmatched==0);
 size_t ax[7],n;uint8_t structural[7],eligible[7];char states[7];OK(saf_coverage_axes(f.c,ax,structural,eligible,7));OK(saf_mapping_states(f.m,states,7));
 size_t fused=SAF_NONE,struct_n=0,fused_n=0;for(size_t i=0;i<7;++i){if(structural[i])++struct_n;if(states[i]=='F'){++fused_n;fused=ax[i];CHECK(!eligible[i]);}}CHECK(struct_n==1&&fused_n==2);
 size_t members[2];OK(saf_coverage_members(f.c,fused,members,2,&n));CHECK(n==2);
 uint64_t split;OK(saf_coverage_projection(f.c,fused,&split,1,&n,nullptr));CHECK(split==3);
 double lengths[7]={NAN,0.25,1,2,0.5,3,4},primitive[7],projected[5],gene[5];uint8_t pv[7],qv[5],gv[5];
 OK(saf_mapping_apply(f.m,lengths,nullptr,7,{primitive,pv,7},{projected,qv,5},{gene,gv,5}));CHECK(projected[fused]==0.75&&qv[fused]);for(auto i:members)CHECK(!pv[i]&&std::isnan(primitive[i]));
 // Rejecting a different topology proposal leaves the accepted immutable plan valid.
 gt[3]=2;gt[5]=1;saf_mapping*proposal=nullptr;OK(saf_mapping_create(f.c,gp,gt,7,&proposal));
 OK(saf_mapping_get_info(proposal,&mi));CHECK(mi.n_reference_absent==1&&mi.n_gene_unmatched==1);
 size_t absent,unmatched;OK(saf_mapping_reference_absent(proposal,&absent,1,&n));CHECK(absent==fused);OK(saf_mapping_gene_unmatched(proposal,&unmatched,1,&n));OK(saf_mapping_gene_split(proposal,unmatched,&split,1,&n,nullptr));CHECK(split==5);
 OK(saf_mapping_states(proposal,states,7));for(auto i:members)CHECK(states[i]=='F');
 saf_support*s=nullptr;OK(saf_support_create(f.r,"gene","fusion",&s));OK(saf_support_update(s,proposal));saf_support_count counts[7];OK(saf_support_read(s,counts,7,&n));for(auto i:members)CHECK(counts[i].n_fused==1&&counts[i].n_dec==0&&counts[i].n_rec==0);
 saf_support_destroy(s);saf_mapping_destroy(proposal);OK(saf_mapping_apply(f.m,lengths,nullptr,7,skip,{projected,qv,5},skip));CHECK(projected[fused]==0.75);
 // Zero/one taxa truth is valid independently of a mapping.
 for(uint64_t small: {UINT64_C(0),UINT64_C(1)}){saf_coverage*c=nullptr;OK(saf_coverage_create(f.r,&small,1,&c));saf_coverage_info ci;OK(saf_coverage_get_info(c,&ci));CHECK(ci.n_projections==0);OK(saf_coverage_axes(c,ax,structural,eligible,7));for(size_t i=0;i<7;++i)CHECK(structural[i]&&ax[i]==SAF_NONE&&!eligible[i]);saf_coverage_destroy(c);}
 // A star free tree has absent reference internal edges; compatibility is not inferred.
 uint64_t full=31;saf_coverage*fullc=nullptr;OK(saf_coverage_create(f.r,&full,1,&fullc));int64_t sp[]={-1,0,0,0,0,0},st[]={-1,0,1,2,3,4};saf_mapping*star=nullptr;OK(saf_mapping_create(fullc,sp,st,6,&star));OK(saf_mapping_get_info(star,&mi));CHECK(mi.n_reference_absent==2&&mi.n_gene_unmatched==0);saf_mapping_destroy(star);saf_coverage_destroy(fullc);
}
void invalid_inputs_and_words(){
 Fixture f;two(f);saf_reference*r=reinterpret_cast<saf_reference*>(1);int64_t p[]={-1,0,0},t[]={-1,0,0};
 int64_t valid_taxon[]={-1,0,1};fail_next_allocation=true;CHECK(saf_reference_create(2,"oom",p,valid_taxon,3,&r)==SAF_NO_MEMORY&&r==nullptr);
 CHECK(saf_reference_create(2,"dup",p,t,3,&r)==SAF_INVALID_ARGUMENT&&r==nullptr);
 t[2]=2;CHECK(saf_reference_create(2,"extra",p,t,3,&r)==SAF_INVALID_ARGUMENT&&r==nullptr);
 t[2]=1;p[2]=9;CHECK(saf_reference_create(2,"range",p,t,3,&r)==SAF_INVALID_ARGUMENT);
 p[2]=-1;CHECK(saf_reference_create(2,"roots",p,t,3,&r)==SAF_INVALID_ARGUMENT);
 int64_t cyc_p[]={-1,0,0,4,3},cyc_t[]={-1,0,1,-1,-1};CHECK(saf_reference_create(2,"cycle",cyc_p,cyc_t,5,&r)==SAF_INVALID_ARGUMENT);
 int64_t unary_p[]={-1,0,1,0},unary_t[]={-1,-1,0,1};CHECK(saf_reference_create(2,"unary",unary_p,unary_t,4,&r)==SAF_INVALID_ARGUMENT);
 uint64_t bad=7;saf_coverage*c=reinterpret_cast<saf_coverage*>(1);CHECK(saf_coverage_create(f.r,&bad,1,&c)==SAF_INVALID_ARGUMENT&&c==nullptr);bad=3;CHECK(saf_coverage_create(f.r,&bad,2,&c)==SAF_INVALID_ARGUMENT);
 bad=1;OK(saf_coverage_create(f.r,&bad,1,&c));p[2]=0;saf_mapping*m=reinterpret_cast<saf_mapping*>(1);CHECK(saf_mapping_create(c,p,t,3,&m)==SAF_INVALID_ARGUMENT&&m==nullptr);saf_coverage_destroy(c);
 for(size_t nt:{63,64,65,127,128,129,302}){
  std::vector<int64_t> pp(nt+1,0),tt(nt+1);pp[0]=-1;tt[0]=-1;for(size_t i=0;i<nt;++i)tt[i+1]=static_cast<int64_t>(i);
  saf_reference*rr=nullptr;OK(saf_reference_create(nt,"boundary-star",pp.data(),tt.data(),pp.size(),&rr));saf_reference_info ri;OK(saf_reference_get_info(rr,&ri));CHECK(ri.n_axes==nt&&ri.n_words==(nt+63)/64);
  std::vector<uint64_t> words(ri.n_words);size_t count;bool seen_high=false;for(size_t a=0;a<ri.n_axes;++a){saf_split_info info;OK(saf_reference_axis(rr,a,words.data(),words.size(),&count,&info));CHECK(info.terminal&&std::min(info.side_count,info.other_count)==1);if(info.side_count==nt-1)seen_high=true;}CHECK(seen_high);
  words.assign(ri.n_words,0);words[0]=1;words[(nt-1)/64]|=UINT64_C(1)<<((nt-1)%64);saf_coverage*cc=nullptr;OK(saf_coverage_create(rr,words.data(),words.size(),&cc));saf_coverage_info ci;OK(saf_coverage_get_info(cc,&ci));CHECK(ci.n_observed==2&&ci.n_projections==1);saf_coverage_destroy(cc);
  if(nt%64){words.back()|=UINT64_C(1)<<(nt%64);CHECK(saf_coverage_create(rr,words.data(),words.size(),&cc)==SAF_INVALID_ARGUMENT);}
  saf_reference_destroy(rr);
 }
}
void support_overflow_and_identity(){
 Fixture f;two(f);saf_support*sum=nullptr,*power=nullptr;OK(saf_support_create(f.r,"gene","overflow",&sum));OK(saf_support_create(f.r,"gene","overflow",&power));OK(saf_support_update(power,f.m));
 for(unsigned i=0;i<64;++i){OK(saf_support_merge(sum,power));if(i<63){saf_support*copy=nullptr;OK(saf_support_create(f.r,"gene","overflow",&copy));OK(saf_support_merge(copy,power));OK(saf_support_merge(power,copy));saf_support_destroy(copy);}}
 size_t count;saf_support_count before{},after{};OK(saf_support_read(sum,&before,1,&count));CHECK(before.n_valid_topologies==UINT64_MAX&&before.n_dec==UINT64_MAX&&before.n_rec==UINT64_MAX);
 CHECK(saf_support_update(sum,f.m)==SAF_COUNT_OVERFLOW);OK(saf_support_read(sum,&after,1,&count));CHECK(std::memcmp(&before,&after,sizeof before)==0);
 CHECK(saf_support_merge(sum,power)==SAF_COUNT_OVERFLOW);OK(saf_support_read(sum,&after,1,&count));CHECK(std::memcmp(&before,&after,sizeof before)==0);
 CHECK(saf_support_merge(sum,sum)==SAF_INVALID_ARGUMENT);saf_support*wrong=nullptr;OK(saf_support_create(f.r,"sample","overflow",&wrong));CHECK(saf_support_merge(sum,wrong)==SAF_INVALID_ARGUMENT);saf_support_destroy(wrong);
 Fixture distinct;two(distinct);OK(saf_support_create(distinct.r,"gene","overflow",&wrong));CHECK(saf_support_update(wrong,f.m)==SAF_INVALID_ARGUMENT);CHECK(saf_support_merge(sum,wrong)==SAF_INVALID_ARGUMENT);saf_support_destroy(wrong);
 saf_support_destroy(sum);saf_support_destroy(power);
}
void thread_errors(){
 std::atomic<int>stage{0};std::string a,b;
 std::thread one([&]{saf_reference_info info;CHECK(saf_reference_get_info(nullptr,&info)==SAF_INVALID_ARGUMENT);a=saf_last_error();stage.fetch_add(1);while(stage.load()<2)std::this_thread::yield();CHECK(a==saf_last_error());});
 std::thread two_thread([&]{Fixture f;two(f);CHECK(saf_reference_get_info(f.r,nullptr)==SAF_INVALID_ARGUMENT);b=saf_last_error();stage.fetch_add(1);while(stage.load()<2)std::this_thread::yield();CHECK(b==saf_last_error());});one.join();two_thread.join();CHECK(!a.empty()&&!b.empty()&&a!=b);
}
int main(){numeric_and_queries();fusion_and_topology();invalid_inputs_and_words();support_overflow_and_identity();thread_errors();std::puts("native_tests: PASS (API/numeric/lifetime/allocation/threads/overflow/boundaries)");}
