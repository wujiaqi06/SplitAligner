#include "splitaligner/fast_core.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <unordered_map>
namespace saf { namespace {
[[noreturn]] void invalid(const char*m){throw Error(SAF_INVALID_ARGUMENT,m);}
size_t pop(const Words&w){size_t n=0;for(auto x:w)n+=static_cast<size_t>(__builtin_popcountll(x));return n;}
bool less_words(const Words&a,const Words&b){for(size_t i=a.size();i-->0;){if(a[i]!=b[i])return a[i]<b[i];}return false;}
Words canonical(Words side,const Words&mask){Words other(mask.size());for(size_t i=0;i<mask.size();++i)other[i]=mask[i]^side[i];return less_words(other,side)?std::move(other):std::move(side);}
struct Hash {size_t operator()(const Words&w)const noexcept {
#ifdef SAF_TEST_FORCE_HASH_COLLISION
 (void)w; return 0;
#else
 size_t h=1469598103934665603ULL;for(auto v:w){h^=static_cast<size_t>(v);h*=1099511628211ULL;h^=static_cast<size_t>(v>>32);}return h;
#endif
}};
using Index=std::unordered_map<Words,size_t,Hash>;
struct Built {Words mask;std::vector<Split> splits;};
Built build_tree(size_t nt,const int64_t*parent,const int64_t*taxon,size_t nn){
 if(nt<2||nn<3||!parent||!taxon)invalid("tree needs at least two taxa and valid node arrays");
 if(nt>static_cast<size_t>(INT64_MAX)||nn>static_cast<size_t>(INT64_MAX))invalid("tree index range exceeds int64");
 size_t width=nt/64+(nt%64!=0); if(width && nn>std::numeric_limits<size_t>::max()/width)invalid("tree bitset size overflow");
 std::vector<std::vector<size_t>> children(nn);size_t root=SAF_NONE;
 for(size_t i=0;i<nn;++i){if(parent[i]==-1){if(root!=SAF_NONE)invalid("multiple roots");root=i;}else{if(parent[i]<0||static_cast<uint64_t>(parent[i])>=nn||static_cast<size_t>(parent[i])==i)invalid("invalid parent index");children[static_cast<size_t>(parent[i])].push_back(i);}}
 if(root==SAF_NONE)invalid("missing root");
 std::vector<size_t> order;order.reserve(nn);order.push_back(root);
 for(size_t k=0;k<order.size();++k){auto u=order[k];for(auto v:children[u])order.push_back(v);}
 if(order.size()!=nn)invalid("tree is disconnected or cyclic");
 Built b;b.mask.assign(width,0);size_t leaves=0;
 for(size_t i=0;i<nn;++i){if(children[i].empty()){if(taxon[i]<0||static_cast<uint64_t>(taxon[i])>=nt)invalid("invalid leaf taxon");auto t=static_cast<size_t>(taxon[i]);uint64_t bit=UINT64_C(1)<<(t%64);if(b.mask[t/64]&bit)invalid("duplicate taxon");b.mask[t/64]|=bit;++leaves;}else {if(taxon[i]!=-1)invalid("internal node carries taxon");if(children[i].size()<2)invalid("unary internal node");}}
 if(leaves<2)invalid("mapping tree needs at least two observed taxa");
 std::vector<Words> sides(nn,Words(width,0));Index index; index.reserve(nn);
 for(auto it=order.rbegin();it!=order.rend();++it){size_t u=*it;
  if(children[u].empty()){size_t t=static_cast<size_t>(taxon[u]);sides[u][t/64]|=UINT64_C(1)<<(t%64);}else for(auto v:children[u])for(size_t j=0;j<width;++j)sides[u][j]|=sides[v][j];
  if(u==root)continue;
  Words side=canonical(sides[u],b.mask);auto hit=index.find(side);
  if(hit==index.end()){Split s;s.side=std::move(side);s.side_count=pop(s.side);s.other_count=leaves-s.side_count;s.terminal=std::min(s.side_count,s.other_count)==1;s.sources.push_back(u);index.emplace(s.side,b.splits.size());b.splits.push_back(std::move(s));}
  else {auto &sources=b.splits[hit->second].sources;if(children[root].size()!=2||parent[u]!=static_cast<int64_t>(root)||sources.size()!=1||parent[sources[0]]!=static_cast<int64_t>(root))invalid("duplicate split is not a binary-root alias");sources.push_back(u);}
 }
 std::sort(b.splits.begin(),b.splits.end(),[](const Split&a,const Split&c){return less_words(a.side,c.side);});
 for(auto&s:b.splits)std::sort(s.sources.begin(),s.sources.end());
 return b;
}
void check_buffer(saf_value_buffer b,size_t needed){
 if(!b.values&&!b.valid){if(b.capacity)invalid("skipped output capacity must be zero");return;}
 if(!b.values||!b.valid)invalid("requested output needs values and valid pointers");
 if(b.capacity<needed)throw Error(SAF_CAPACITY,"apply output capacity too small");
}
struct Range{uintptr_t begin,end;};
Range range(const void*p,size_t n,size_t width){if(!p||!n)return {0,0};if(n>std::numeric_limits<size_t>::max()/width)invalid("buffer byte size overflow");auto start=reinterpret_cast<uintptr_t>(p);auto bytes=n*width;if(start>std::numeric_limits<uintptr_t>::max()-bytes)invalid("buffer address range overflow");return {start,start+bytes};}
bool overlaps(Range a,Range b){return a.begin<a.end&&b.begin<b.end&&a.begin<b.end&&b.begin<a.end;}
void aggregate(const Split*s,const double*values,const uint8_t*valid,double&out,uint8_t&available){
 out=std::numeric_limits<double>::quiet_NaN();available=0;if(!s)return;
 for(size_t source:s->sources){if(valid&&!valid[source])return;if(!std::isfinite(values[source]))return;}
 double total=0;for(size_t source:s->sources)total+=values[source];
 if(!std::isfinite(total))throw Error(SAF_NUMERIC_OVERFLOW,"finite source aggregation overflow");
 out=total;available=1;
}
bool add_ok(uint64_t a,uint64_t b){return b<=UINT64_MAX-a;}
bool count_ok(const saf_support_count&a,const saf_support_count&b){return add_ok(a.n_valid_topologies,b.n_valid_topologies)&&add_ok(a.n_struct,b.n_struct)&&add_ok(a.n_fused,b.n_fused)&&add_ok(a.n_dec,b.n_dec)&&add_ok(a.n_rec,b.n_rec)&&add_ok(a.n_not_recovered,b.n_not_recovered);}
void add_count(saf_support_count&a,const saf_support_count&b){a.n_valid_topologies+=b.n_valid_topologies;a.n_struct+=b.n_struct;a.n_fused+=b.n_fused;a.n_dec+=b.n_dec;a.n_rec+=b.n_rec;a.n_not_recovered+=b.n_not_recovered;}
saf_support_count increment(char state,bool eligible){saf_support_count x{};x.n_valid_topologies=1;x.n_struct=state=='S';x.n_fused=state=='F';x.n_dec=eligible;x.n_rec=eligible&&state=='D';x.n_not_recovered=eligible&&state=='U';return x;}
}
std::shared_ptr<const Reference> Reference::create(size_t nt,const char*universe,const int64_t*p,const int64_t*t,size_t nn){
 if(!universe||!*universe)invalid("universe_id must be nonempty");
 Built b=build_tree(nt,p,t,nn);
 if(pop(b.mask)!=nt)invalid("reference must contain every taxon");
 std::shared_ptr<Reference>r(new Reference);r->n_taxa_=nt;r->n_nodes_=nn;r->full_=std::move(b.mask);r->axes_=std::move(b.splits);r->universe_=universe;return r;
}
std::shared_ptr<const Coverage> Coverage::create(std::shared_ptr<const Reference>r,const uint64_t*mask,size_t nw){
 if(!r||nw!=r->n_words()||!mask)invalid("coverage width or reference invalid");
 for(size_t i=0;i<nw;++i)if(mask[i]&~r->full()[i])invalid("coverage has nonzero padding or extra taxa");
 std::shared_ptr<Coverage>c(new Coverage);c->reference_=std::move(r);c->mask_.assign(mask,mask+nw);c->observed_=pop(c->mask_);size_t n=c->reference_->axes().size();c->axis_to_projection_.assign(n,SAF_NONE);c->structural_.assign(n,0);c->eligible_.assign(n,0);Index index;index.reserve(n);
 for(size_t a=0;a<n;++a){Words side(nw);for(size_t j=0;j<nw;++j)side[j]=c->reference_->axes()[a].side[j]&mask[j];size_t cnt=pop(side);if(cnt==0||cnt==c->observed_){c->structural_[a]=1;continue;}side=canonical(std::move(side),c->mask_);auto hit=index.find(side);size_t id;
  if(hit==index.end()){id=c->projections_.size();Projection p;p.split.side=std::move(side);p.split.side_count=pop(p.split.side);p.split.other_count=c->observed_-p.split.side_count;p.split.terminal=std::min(p.split.side_count,p.split.other_count)==1;index.emplace(p.split.side,id);c->projections_.push_back(std::move(p));}else id=hit->second;c->projections_[id].members.push_back(a);
 }
 std::sort(c->projections_.begin(),c->projections_.end(),[](const Projection&a,const Projection&b){return less_words(a.split.side,b.split.side);});
 for(size_t i=0;i<c->projections_.size();++i){const auto&p=c->projections_[i];for(size_t a:p.members){c->axis_to_projection_[a]=i;c->eligible_[a]=(p.members.size()==1&&(c->reference_->axes()[a].terminal||!p.split.terminal));}}
 return c;
}
std::shared_ptr<const Mapping> Mapping::create(std::shared_ptr<const Coverage>c,const int64_t*p,const int64_t*t,size_t nn){
 if(!c)invalid("null coverage");
 Built b=build_tree(c->reference()->n_taxa(),p,t,nn);
 if(b.mask!=c->mask())invalid("gene leafset does not equal coverage");
 std::shared_ptr<Mapping>m(new Mapping);m->coverage_=std::move(c);m->n_nodes_=nn;m->gene_splits_=std::move(b.splits);const auto&projs=m->coverage_->projections();m->projection_to_gene_.assign(projs.size(),SAF_NONE);m->gene_to_projection_.assign(m->gene_splits_.size(),SAF_NONE);m->recovered_.assign(projs.size(),0);Index index;index.reserve(m->gene_splits_.size());for(size_t i=0;i<m->gene_splits_.size();++i)index.emplace(m->gene_splits_[i].side,i);
 for(size_t i=0;i<projs.size();++i){auto hit=index.find(projs[i].split.side);if(hit==index.end())m->reference_absent_.push_back(i);else{m->projection_to_gene_[i]=hit->second;m->gene_to_projection_[hit->second]=i;m->recovered_[i]=1;}}
 for(size_t i=0;i<m->gene_splits_.size();++i)if(m->gene_to_projection_[i]==SAF_NONE)m->gene_unmatched_.push_back(i);
 m->states_.resize(m->coverage_->axis_to_projection().size());for(size_t a=0;a<m->states_.size();++a){size_t q=m->coverage_->axis_to_projection()[a];char state='S';if(q!=SAF_NONE){if(projs[q].members.size()>1)state='F';else state=m->recovered_[q]?'D':'U';}m->states_[a]=state;}
 return m;
}
void Mapping::apply(const double*values,const uint8_t*valid,size_t nn,saf_value_buffer primitive,saf_value_buffer projection,saf_value_buffer gene)const{
 const size_t sizes[3]={states_.size(),projection_to_gene_.size(),gene_splits_.size()};saf_value_buffer out[3]={primitive,projection,gene};bool requested=false;for(size_t i=0;i<3;++i){check_buffer(out[i],sizes[i]);requested|=out[i].values!=nullptr;}if(!requested)return;
 if(nn!=n_nodes_||!values)invalid("apply needs original node-layout lengths and exact node count");
 Range in[2]={range(values,nn,sizeof(double)),range(valid,nn,sizeof(uint8_t))};Range output[6];for(size_t i=0;i<3;++i){output[2*i]=range(out[i].values,sizes[i],sizeof(double));output[2*i+1]=range(out[i].valid,sizes[i],sizeof(uint8_t));}for(size_t i=0;i<6;++i){for(auto r:in)if(overlaps(output[i],r))invalid("input/output buffers overlap");for(size_t j=0;j<i;++j)if(overlaps(output[i],output[j]))invalid("output buffers overlap");}
 if(primitive.values)for(size_t a=0;a<sizes[0];++a){const Split*s=nullptr;if(states_[a]=='D'&&coverage_->eligible()[a])s=&gene_splits_[projection_to_gene_[coverage_->axis_to_projection()[a]]];aggregate(s,values,valid,primitive.values[a],primitive.valid[a]);}
 if(projection.values)for(size_t i=0;i<sizes[1];++i){size_t edge=projection_to_gene_[i];aggregate(edge==SAF_NONE?nullptr:&gene_splits_[edge],values,valid,projection.values[i],projection.valid[i]);}
 if(gene.values)for(size_t i=0;i<sizes[2];++i)aggregate(&gene_splits_[i],values,valid,gene.values[i],gene.valid[i]);
}
std::unique_ptr<Support> Support::create(std::shared_ptr<const Reference>r,const char*unit,const char*cohort){if(!r||!unit||!*unit||!cohort||!*cohort)invalid("support reference, unit and cohort required");std::unique_ptr<Support>s(new Support);s->reference_=std::move(r);s->unit_=unit;s->cohort_=cohort;s->counts_.resize(s->reference_->axes().size());return s;}
void Support::update(const Mapping&m){if(m.coverage()->reference().get()!=reference_.get())invalid("support reference identity mismatch");for(size_t a=0;a<counts_.size();++a)if(!count_ok(counts_[a],increment(m.states()[a],m.coverage()->eligible()[a])))throw Error(SAF_COUNT_OVERFLOW,"support update counter overflow");for(size_t a=0;a<counts_.size();++a)add_count(counts_[a],increment(m.states()[a],m.coverage()->eligible()[a]));}
void Support::merge(const Support&s){if(this==&s)invalid("support self merge is not a disjoint shard");if(reference_.get()!=s.reference_.get()||unit_!=s.unit_||cohort_!=s.cohort_)invalid("support merge identity, unit or cohort mismatch");for(size_t a=0;a<counts_.size();++a)if(!count_ok(counts_[a],s.counts_[a]))throw Error(SAF_COUNT_OVERFLOW,"support merge counter overflow");for(size_t a=0;a<counts_.size();++a)add_count(counts_[a],s.counts_[a]);}
}
