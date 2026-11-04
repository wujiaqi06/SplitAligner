#ifndef SPLITALIGNER_FAST_CORE_HPP
#define SPLITALIGNER_FAST_CORE_HPP
#include "fast_c_api.h"
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace saf {
using Words=std::vector<uint64_t>;
struct Error:std::runtime_error { saf_status code; Error(saf_status c,const char *m):std::runtime_error(m),code(c){} };
struct Split { Words side; size_t side_count=0, other_count=0; bool terminal=false; std::vector<size_t> sources; };
struct Projection { Split split; std::vector<size_t> members; };
class Reference final {
public:
 Reference(const Reference&)=delete;
 Reference(Reference&&)=delete;
 Reference& operator=(const Reference&)=delete;
 Reference& operator=(Reference&&)=delete;
 static std::shared_ptr<const Reference> create(size_t,const char*,const int64_t*,const int64_t*,size_t);
 size_t n_taxa()const{return n_taxa_;} size_t n_nodes()const{return n_nodes_;} size_t n_words()const{return full_.size();}
 const Words& full()const{return full_;} const std::vector<Split>& axes()const{return axes_;} const std::string& universe()const{return universe_;}
private:
 Reference()=default;
 size_t n_taxa_=0,n_nodes_=0; Words full_; std::vector<Split> axes_; std::string universe_;
};
class Coverage final {
public:
 Coverage(const Coverage&)=delete;
 Coverage(Coverage&&)=delete;
 Coverage& operator=(const Coverage&)=delete;
 Coverage& operator=(Coverage&&)=delete;
 static std::shared_ptr<const Coverage> create(std::shared_ptr<const Reference>,const uint64_t*,size_t);
 const std::shared_ptr<const Reference>& reference()const{return reference_;}
 const Words& mask()const{return mask_;} size_t observed()const{return observed_;}
 const std::vector<Projection>& projections()const{return projections_;}
 const std::vector<size_t>& axis_to_projection()const{return axis_to_projection_;}
 const std::vector<uint8_t>& structural()const{return structural_;}
 const std::vector<uint8_t>& eligible()const{return eligible_;}
private:
 Coverage()=default;
 std::shared_ptr<const Reference> reference_; Words mask_; size_t observed_=0;
 std::vector<Projection> projections_; std::vector<size_t> axis_to_projection_; std::vector<uint8_t> structural_,eligible_;
};
class Mapping final {
public:
 Mapping(const Mapping&)=delete;
 Mapping(Mapping&&)=delete;
 Mapping& operator=(const Mapping&)=delete;
 Mapping& operator=(Mapping&&)=delete;
 static std::shared_ptr<const Mapping> create(std::shared_ptr<const Coverage>,const int64_t*,const int64_t*,size_t);
 const std::shared_ptr<const Coverage>& coverage()const{return coverage_;} size_t n_nodes()const{return n_nodes_;}
 const std::vector<Split>& gene_splits()const{return gene_splits_;}
 const std::vector<size_t>& projection_to_gene()const{return projection_to_gene_;}
 const std::vector<size_t>& gene_to_projection()const{return gene_to_projection_;}
 const std::vector<size_t>& reference_absent()const{return reference_absent_;}
 const std::vector<size_t>& gene_unmatched()const{return gene_unmatched_;}
 const std::vector<char>& states()const{return states_;}
 const std::vector<uint8_t>& recovered()const{return recovered_;}
 void apply(const double*,const uint8_t*,size_t,saf_value_buffer,saf_value_buffer,saf_value_buffer)const;
private:
 Mapping()=default;
 std::shared_ptr<const Coverage> coverage_; size_t n_nodes_=0; std::vector<Split> gene_splits_;
 std::vector<size_t> projection_to_gene_,gene_to_projection_,reference_absent_,gene_unmatched_; std::vector<char> states_; std::vector<uint8_t> recovered_;
};
class Support final {
public:
 static std::unique_ptr<Support> create(std::shared_ptr<const Reference>,const char*,const char*);
 void update(const Mapping&); void merge(const Support&);
 const std::vector<saf_support_count>& counts()const{return counts_;}
private:
 Support()=default;
 std::shared_ptr<const Reference> reference_; std::string unit_,cohort_; std::vector<saf_support_count> counts_;
};
}
#endif
