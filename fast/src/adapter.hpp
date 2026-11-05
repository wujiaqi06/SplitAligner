#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <locale.h>
#include <map>
#include <limits>
#include <memory>
#include <openssl/evp.h>
#include <openssl/sha.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace safio {
inline bool space(char c){return c==' '||c=='\t'||c=='\r'||c=='\n'||c=='\f'||c=='\v';}
inline bool blank(const std::string&s){return std::all_of(s.begin(),s.end(),space);}
// Return the size of one Unicode scalar in UTF-8, or zero for malformed input.
inline size_t utf8_width(const std::string& s, size_t i) {
 const unsigned char c=s[i];
 if(c<0x80)return 1;
 size_t width;
 if(c>=0xc2 && c<=0xdf)width=2;
 else if(c>=0xe0 && c<=0xef)width=3;
 else if(c>=0xf0 && c<=0xf4)width=4;
 else return 0;
 if(s.size()-i<width)return 0;
 for(size_t k=1;k<width;++k)if((static_cast<unsigned char>(s[i+k])&0xc0)!=0x80)return 0;
 const unsigned char second=s[i+1];
 if((c==0xe0 && second<0xa0)||(c==0xed && second>=0xa0)||
    (c==0xf0 && second<0x90)||(c==0xf4 && second>=0x90))return 0;
 return width;
}
inline bool utf8(const std::string& s) {
 for(size_t i=0;i<s.size();){auto n=utf8_width(s,i);if(!n)return false;i+=n;}
 return true;
}
inline void strip_initial_bom(std::string& s) {
 if(s.compare(0,3,"\xef\xbb\xbf")==0)s.erase(0,3);
}
inline void identity(const std::string&s,const char*w){if(s.empty()||blank(s)||!utf8(s))throw std::runtime_error(std::string("invalid ")+w);for(unsigned char c:s)if(c<32||c==127)throw std::runtime_error(std::string("control character in ")+w);}
inline std::string json(const std::string& s) {
 std::string out="\"";
 for(size_t i=0;i<s.size();) {
  const unsigned char c=s[i];
  if(c>=0x80) {
   const size_t n=utf8_width(s,i);
   if(!n){out+="\\ufffd";++i;continue;}
   if(n==3 && s.compare(i,3,"\xe2\x80\xa8")==0)out+="\\u2028";
   else if(n==3 && s.compare(i,3,"\xe2\x80\xa9")==0)out+="\\u2029";
   else out.append(s,i,n);
   i+=n;continue;
  }
  ++i;
  switch(c){
   case '"':out+="\\\"";break;case '\\':out+="\\\\";break;
   case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;
   default:if(c<32){const char* h="0123456789abcdef";out+="\\u00";out+=h[c>>4];out+=h[c&15];}
           else out+=char(c);
  }
 }
 return out+'"';
}
inline std::string hex_bytes(const unsigned char*p,size_t n){const char*h="0123456789abcdef";std::string s;for(size_t i=0;i<n;++i){s+=h[p[i]>>4];s+=h[p[i]&15];}return s;}
inline std::string digest(const std::string&s){unsigned char h[32];SHA256((const unsigned char*)s.data(),s.size(),h);return hex_bytes(h,32);}
class InputDigest {
 std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> ctx{EVP_MD_CTX_new(), EVP_MD_CTX_free};
 uintmax_t size_=0;
public:
 InputDigest() {
  if(!ctx || EVP_DigestInit_ex(ctx.get(),EVP_sha256(),nullptr)!=1)throw std::runtime_error("SHA-256 initialization failed");
 }
 void add(const char* data,size_t n) {
  if(n>std::numeric_limits<uintmax_t>::max()-size_)throw std::runtime_error("input byte count overflow");
  if(EVP_DigestUpdate(ctx.get(),data,n)!=1)throw std::runtime_error("SHA-256 update failed");
  size_+=n;
 }
 uintmax_t size()const{return size_;}
 std::string finish() {
  unsigned char hash[EVP_MAX_MD_SIZE];unsigned int n=0;
  if(EVP_DigestFinal_ex(ctx.get(),hash,&n)!=1 || n!=32)throw std::runtime_error("SHA-256 finalization failed");
  return hex_bytes(hash,n);
 }
};
template<class Consumer>inline void read_chunks(const std::string& path, Consumer consume) {
 std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("open: "+path);
 char buffer[65536];
 while(f.read(buffer,sizeof buffer)||f.gcount())consume(buffer,static_cast<size_t>(f.gcount()));
 if(f.bad()||!f.eof())throw std::runtime_error("read: "+path);
}
inline std::string file_digest(const std::string& path) {
 InputDigest hash;read_chunks(path,[&](const char* b,size_t n){hash.add(b,n);});return hash.finish();
}
inline std::string bits_hex(const std::vector<uint64_t>&w){std::ostringstream s;s<<std::hex<<std::setfill('0');for(auto i=w.rbegin();i!=w.rend();++i)s<<std::setw(16)<<*i;return s.str();}
inline bool missing(const std::string&s){std::string u;for(unsigned char c:s)u+=c>='a'&&c<='z'?char(c-32):char(c);if(u=="NA"||u=="N/A"||u=="."||u=="?"||u=="NULL"||u=="NONE")return true;if(!u.empty()&&(u[0]=='+'||u[0]=='-'))u.erase(u.begin());return u=="NAN"||u=="INF"||u=="INFINITY";}
inline bool decimal(const std::string&s){size_t i=0;if(i<s.size()&&(s[i]=='+'||s[i]=='-'))++i;size_t d=0;while(i<s.size()&&s[i]>='0'&&s[i]<='9'){++i;++d;}if(i<s.size()&&s[i]=='.'){++i;while(i<s.size()&&s[i]>='0'&&s[i]<='9'){++i;++d;}}if(!d)return false;if(i<s.size()&&(s[i]=='e'||s[i]=='E')){++i;if(i<s.size()&&(s[i]=='+'||s[i]=='-'))++i;size_t start=i;while(i<s.size()&&s[i]>='0'&&s[i]<='9')++i;if(i==start)return false;}return i==s.size();}
inline double number(const std::string&s){struct Locale{locale_t p=newlocale(LC_NUMERIC_MASK,"C",nullptr);~Locale(){if(p)freelocale(p);}};static Locale loc;if(!loc.p)throw std::runtime_error("C locale unavailable");char*end=nullptr;double v=strtod_l(s.c_str(),&end,loc.p);if(end!=s.data()+s.size()||!std::isfinite(v))throw std::runtime_error("numeric range: "+s);bool nz=false;for(char c:s){if(c=='e'||c=='E')break;if(c>='1'&&c<='9')nz=true;}if(v==0&&nz)throw std::runtime_error("numeric underflow to zero: "+s);return v;}
struct Node{int64_t parent=-1;std::vector<size_t> children;std::string label,length;std::vector<std::string> comments;bool has_length=false;bool tip_was_quoted=false;};
struct Tree{std::vector<Node> nodes;size_t root=0;};
class Parser {
 const std::string&s;size_t p=0;Tree t;
 void skip(std::vector<std::string>&comments){while(p<s.size()){if(space(s[p])){++p;continue;}if(s[p]!='[')break;size_t start=++p;int d=1;while(p<s.size()&&d){if(s[p]=='[')++d;else if(s[p]==']')--d;++p;}if(d)throw std::runtime_error("unclosed comment");comments.push_back(s.substr(start,p-start-1));}}
 std::string label(){if(p==s.size())throw std::runtime_error("missing label");std::string v;if(s[p]=='\''){++p;bool done=false;while(p<s.size()){char c=s[p++];if(c=='\''){if(p<s.size()&&s[p]=='\''){v+='\'';++p;}else{done=true;break;}}else v+=c;}if(!done)throw std::runtime_error("unclosed quote");}else{size_t start=p;while(p<s.size()&&!space(s[p])&&std::string("(),:;[]").find(s[p])==std::string::npos)++p;v=s.substr(start,p-start);}identity(v,"taxon/annotation label");return v;}
 size_t subtree(int64_t parent,size_t depth){
  if(depth>4096)throw std::runtime_error("Newick nesting exceeds supported depth 4096");std::vector<std::string> pre;skip(pre);if(p==s.size())throw std::runtime_error("missing subtree");size_t v=t.nodes.size();t.nodes.emplace_back();t.nodes[v].parent=parent;t.nodes[v].comments=std::move(pre);
  if(s[p]=='('){++p;for(;;){size_t c=subtree((int64_t)v,depth+1);t.nodes[v].children.push_back(c);skip(t.nodes[v].comments);if(p==s.size())throw std::runtime_error("missing right parenthesis");if(s[p]==','){++p;continue;}if(s[p]==')'){++p;break;}throw std::runtime_error("expected comma or closing parenthesis");}if(t.nodes[v].children.size()<2)throw std::runtime_error("unary/empty internal node");skip(t.nodes[v].comments);if(p<s.size()&&std::string(":,);").find(s[p])==std::string::npos)t.nodes[v].label=label();}
  else{if(std::string(",):;]").find(s[p])!=std::string::npos)throw std::runtime_error("empty or invalid subtree");t.nodes[v].tip_was_quoted=s[p]=='\'';t.nodes[v].label=label();}
  skip(t.nodes[v].comments);if(p<s.size()&&s[p]==':'){++p;skip(t.nodes[v].comments);size_t start=p;while(p<s.size()&&!space(s[p])&&std::string("(),:;[]").find(s[p])==std::string::npos)++p;auto token=s.substr(start,p-start);if(!missing(token)&&!decimal(token))throw std::runtime_error("invalid length token: "+token);t.nodes[v].length=token;t.nodes[v].has_length=true;skip(t.nodes[v].comments);}
  if(p<s.size()&&std::string(",);").find(s[p])==std::string::npos)throw std::runtime_error("invalid token order or duplicate length");return v;
 }
public:explicit Parser(const std::string&text):s(text){}Tree parse(){if(!utf8(s))throw std::runtime_error("input is not valid UTF-8");t.root=subtree(-1,0);if(p==s.size()||s[p++]!=';')throw std::runtime_error("missing Newick semicolon");skip(t.nodes[t.root].comments);if(p!=s.size())throw std::runtime_error("text after tree");return std::move(t);}
};
inline std::vector<std::string> leaves(const Tree&t){std::vector<std::string> n;for(auto&v:t.nodes)if(v.children.empty())n.push_back(v.label);std::sort(n.begin(),n.end(),[](const std::string&a,const std::string&b){return std::lexicographical_compare(a.begin(),a.end(),b.begin(),b.end(),[](char x,char y){return (unsigned char)x<(unsigned char)y;});});if(std::adjacent_find(n.begin(),n.end())!=n.end())throw std::runtime_error("duplicate taxon");return n;}
inline void arrays(const Tree&t,const std::map<std::string,size_t>&ids,std::vector<int64_t>&parent,std::vector<int64_t>&taxon,std::vector<uint64_t>&mask){parent.resize(t.nodes.size());taxon.assign(t.nodes.size(),-1);mask.assign((ids.size()+63)/64,0);for(size_t i=0;i<t.nodes.size();++i){parent[i]=t.nodes[i].parent;if(t.nodes[i].children.empty()){auto it=ids.find(t.nodes[i].label);if(it==ids.end())throw std::runtime_error("unknown taxon: "+t.nodes[i].label);size_t id=it->second;if(mask[id/64]&(uint64_t(1)<<(id%64)))throw std::runtime_error("duplicate taxon");mask[id/64]|=uint64_t(1)<<(id%64);taxon[i]=(int64_t)id;}}}
inline void lengths(const Tree&t,std::vector<double>&values,std::vector<uint8_t>&valid){values.assign(t.nodes.size(),0);valid.assign(t.nodes.size(),0);for(size_t i=0;i<t.nodes.size();++i)if(i!=t.root&&t.nodes[i].has_length&&!missing(t.nodes[i].length)){values[i]=number(t.nodes[i].length);valid[i]=1;}}
inline std::pair<std::string,std::string> gene_line(std::string line,bool prefix){if(!line.empty()&&line.back()=='\r')line.pop_back();size_t p=line.find(prefix?'(':'\t');if(p==std::string::npos)throw std::runtime_error("gene line does not match explicit gene format");std::string id=line.substr(0,p);identity(id,"gene ID");return {id,line.substr(p+(prefix?0:1))};}
}
