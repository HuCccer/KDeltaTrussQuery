#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<utility>
#include<algorithm>
#include<iomanip> 
#include<queue> 
#include<cstdlib>
#include<numeric>
#include<unordered_map>
#include<limits.h>
#include<string.h>
#include<stdint.h>
#include<random>

using namespace std;

// edge type
typedef std::pair<uint32_t, uint32_t> EdgT;
// triangle type
typedef std::tuple<uint32_t, uint32_t, uint32_t> Triangle;
typedef vector<uint32_t> VI;

struct ArrayEntry{
    uint32_t vid;
    uint32_t eid;
};

vector<unordered_map<uint32_t, bool>> Mc;

// data members
uint32_t n_;  // the # of vertices
uint32_t m_;  // the # of edges
uint32_t k_max;	//the # maximum k value
uint32_t t_max; //the #maximum timestamps
// the adjacency array representation
vector<vector<ArrayEntry>> adj_;
// the truss numbers
VI k_;
// the k-support
VI ks_;
//the support
VI s_;
// the edge peeling order
VI ord_;
//the position of edges in order 
VI pos;
// the set of edges
vector<EdgT> edges_;
//the timestamps of edges
vector<VI> edge2times;
//delta-triangle list
vector<vector<Triangle>> time2triangle(4000);
//TC-index
vector<vector<pair<uint32_t,uint32_t>>> TC_index;
vector<vector<uint32_t>> TC;
//record the length of TC-index
vector<int> kTClen;

time_t t_start, t_end;
uint32_t mine , maxe;

class TemporalGraph {
    public:
    uint32_t timestamps;
    uint32_t vertices;
    uint32_t temporal_edges;

    TemporalGraph(string dataset_path) {
        ifstream in(dataset_path);
        in >> timestamps >> vertices >> temporal_edges;
		print_info(dataset_path);
		cout << endl << "Load Graph.." << endl;
		uint32_t _t, v1, v2;
		n_ = vertices;
		adj_.resize(n_);
		map<EdgT,vector<uint32_t>> ET;
        for(uint32_t i = 0; i < temporal_edges; ++i) {   
            in >> _t >> v1 >> v2;
			ET[{v1,v2}].emplace_back(_t);
		}
        in.close();
		m_ = ET.size();
		edges_.resize(m_);
		edge2times.resize(m_);
		Mc.resize(m_);
		uint32_t edgeIndex = 0;
		for(auto it=ET.begin(); it!=ET.end(); ++it){
			adj_[it->first.first].emplace_back(ArrayEntry{it->first.second,edgeIndex});
			adj_[it->first.second].emplace_back(ArrayEntry{it->first.first,edgeIndex});
			edges_[edgeIndex] = it->first;
			edge2times[edgeIndex] = it->second;
			edgeIndex++;
		}
		for (uint32_t vid = 0; vid < n_; ++vid) adj_[vid].shrink_to_fit();
		// clear ET
  		decltype(ET)().swap(ET);
		cout << "Load graph successfully!" << endl;
    }

    private:
		void print_info(string dataset_path) {
			cout << "-----------info--------------" << endl;
			cout << "Number of timestamps: " << timestamps << endl;
			cout << "Number of vertices: " << vertices << endl;
			cout << "Number of temporal edges: " << temporal_edges << endl;
			cout << "-----------------------------" << endl;
		}
};

inline bool Pred(uint32_t i, uint32_t j) {
    return adj_[i].size() > adj_[j].size() || (adj_[i].size() == adj_[j].size() && i > j);
}

inline void removeEdge(uint32_t eid){
	uint32_t v1 = edges_[eid].first;
    uint32_t v2 = edges_[eid].second;
    size_t p1 = 0, p2 = 0;
	while (adj_[v1][p1].vid != v2) ++p1;
    while (adj_[v2][p2].vid != v1) ++p2;
	adj_[v1].erase(adj_[v1].begin() + p1);
    adj_[v2].erase(adj_[v2].begin() + p2);
}

inline void getminmax(uint32_t a, uint32_t b, uint32_t c){
	if (a > b) a ^= b ^= a ^= b;
	if (a > c) a ^= c ^= a ^= c;
	if (b > c){maxe = b; mine = a;}
	else {maxe = c; mine = a;}
}

inline Triangle sortTri(uint32_t a, uint32_t b, uint32_t c){
	if (a > b) a ^= b ^= a ^= b;
	if (a > c) a ^= c ^= a ^= c;
	if (b > c) b ^= c ^= b ^= c;
	return {a, b, c};
}

uint32_t get_mst(const VI &times1,const VI &times2,const VI &times3){
	uint32_t i = 0, j = 0, k = 0;
	uint32_t n = times1.size(), m = times2.size(), p = times3.size();
	uint32_t delta_min = 20000;
	while(i < n && j < m && k < p){
		uint32_t minvalue = min({(uint32_t)times1[i], times2[j], times3[k]});
		uint32_t delta = max({(uint32_t)times1[i], times2[j], times3[k]}) - minvalue;
		if(delta < delta_min) delta_min = delta;
		if(delta_min == 0) return delta_min;
		if(minvalue == times1[i]) i++;
		else if(minvalue == times2[j]) j++;
		else k++;
	}
	return delta_min;
}

void print_exec_time() {
	cout << "-----------------------------" << endl; 
	double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC); 
    cout << "Execution time: " << fixed << time_taken << setprecision(5); 
    cout << " sec " << endl;
	cout << "-----------------------------" << endl; 
}


void print_info(int t_s, int N) {
	int MOD_INFO = 5;
	if(t_s % MOD_INFO == 0) cerr << t_s << " of " << N << endl; 
}

void preProcess(){
	// 1. initialization
  	// 1.1. sort the vertices in non-ascending order of degree
 	VI verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), Pred);
  	// 1.2. call the "forward" algorithm to list triangles
  	k_.resize(m_);s_.resize(m_);
  	vector<vector<ArrayEntry>> A(n_);
	for (uint32_t v : verts) {
    	for (auto ae : adj_[v]) {
			uint32_t u = ae.vid;
			uint32_t e = ae.eid;
			if (!Pred(v, u)) continue;
			size_t pv = 0, pu = 0;
			while (pv < A[v].size() && pu < A[u].size()) {
				if (A[v][pv].vid == A[u][pu].vid) {
					++k_[A[v][pv].eid]; ++k_[A[u][pu].eid];
					++k_[e];
					uint32_t interval = get_mst(edge2times[e],edge2times[A[v][pv].eid],edge2times[A[u][pu].eid]);
					if(interval > t_max) t_max = interval;
					time2triangle[interval].emplace_back(sortTri(e,A[v][pv].eid,A[u][pu].eid));
					++pv; ++pu;
				} else if (Pred(A[v][pv].vid, A[u][pu].vid)) {
				++pv;
				} else {
				++pu;
				}
			}
      		A[u].emplace_back(ArrayEntry{v, e});
		}
	}
	// 2. decomposition
	// 2.1. sort the edges according to their supports
	memcpy(s_.data(),k_.data(),s_.size() * sizeof(uint32_t));
	const uint32_t max_sup = *max_element(k_.cbegin(), k_.cend());
	VI bin(max_sup + 1, 0);
  	for (uint32_t eid = 0; eid < m_; ++eid) ++bin[k_[eid]];
	for (uint32_t i = 0, start = 0; i <= max_sup; ++i) {
		start += bin[i];
		bin[i] = start - bin[i];
	}

  	ord_.resize(m_);
	pos.resize(m_);
  	for (uint32_t eid = 0; eid < m_; ++eid) {
    	pos[eid] = bin[k_[eid]];
    	ord_[pos[eid]] = eid;
    	++bin[k_[eid]];
  	}
	rotate(bin.rbegin(), bin.rbegin() + 1, bin.rend());
  	bin[0] = 0;
  	// 2.2. peeling
  	ks_.resize(m_, 0);
	vector<bool> removed(m_, false);
  	uint32_t k = 0;
  	for (uint32_t i = 0; i < m_; ++i) {
    	k = max(k, k_[ord_[i]]);
    	// ASSERT(bin[k] == i);
    	const uint32_t eid = ord_[i];
    	++bin[k_[eid]];
    	removed[eid] = true;
    	// find triangles containing the edge with ID eid
    	vector<std::pair<uint32_t, uint32_t>> tris; {
      		const uint32_t v1 = edges_[eid].first;
      		const uint32_t v2 = edges_[eid].second;
      		size_t p1 = 0, p2 = 0;
      		while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
				if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
					tris.emplace_back(make_pair(adj_[v1][p1].eid, adj_[v2][p2].eid));
					++p1; ++p2;
				} else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
					++p1;
				} else {
					++p2;
				}
      		}
    	}
		// update ks_[eid]
		for (const auto tri : tris) {
			const uint32_t e1 = tri.first;
			const uint32_t e2 = tri.second;
			if (k_[e1] >= k && k_[e2] >= k) ++ks_[eid];
			if (removed[e1] || removed[e2]) continue;
			for (const uint32_t e : {e1, e2}) {
				if (k_[e] > k) {
				const uint32_t pe3 = bin[k_[e]];
				const uint32_t pe = pos[e];
				if (pe3 != pe) {
					const uint32_t e3 = ord_[pe3];
					ord_[pe] = e3;
					pos[e3] = pe;
					ord_[pe3] = e;
					pos[e] = pe3;
				}
				++bin[k_[e]];
				--k_[e];
				}
			}
		}
  	}
	k_max = k;
	TC.resize(k_max+1);
	TC_index.resize(k_max+1);
	kTClen.resize(k_max+1);
}


void tc_Index_Construction(string dataset_path){
	TemporalGraph tgraph = TemporalGraph(dataset_path);
	t_start = clock();
	preProcess();
	queue<uint32_t> q;
	vector<bool> Ins(m_,false);
	for(uint32_t t = t_max ; t > 0 ; --t){
		if(time2triangle[t].size() == 0) continue;
		for(auto tri:time2triangle[t]){
			uint32_t e1 = std::get<0>(tri), e2 = std::get<1>(tri), e3 = std::get<2>(tri);
			Mc[e1][e3] = true;
			//update support
			s_[e1]--; s_[e2]--; s_[e3]--;
            uint32_t k1 = k_[e1], k2 = k_[e2], k3 = k_[e3];
			//update ks
			if(k1<=k2 && k1<=k3) {
				if(ks_[e1]-- == k1) { q.push(e1); Ins[e1] = true;}
			} 
			if(k2<=k1 && k2<=k3) {
				if(ks_[e2]-- == k2) { q.push(e2); Ins[e2] = true;}
			} 
			if(k3<=k1 && k3<=k2) {
				if(ks_[e3]-- == k3) { q.push(e3); Ins[e3] = true;}
			}
			//update the trussness of edge
			while(!q.empty()){
				uint32_t e=q.front();
				q.pop();
                Ins[e] = false;
				uint32_t k = k_[e];
				uint32_t supe = s_[e];
				if(supe == 0){
					TC[1].emplace_back(e);
                    k_[e] = 0;
					removeEdge(e);
					continue;
				}
				uint32_t nfound=0;
				uint32_t kse = 0;
				uint32_t v1 = edges_[e].first;
      			uint32_t v2 = edges_[e].second;
      			size_t p1 = 0, p2 = 0;
      			while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
					if(nfound == supe) break;
					if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
						getminmax(e, adj_[v1][p1].eid, adj_[v2][p2].eid);
						if(Mc[mine].find(maxe) != Mc[mine].end()) {
							++p1; ++p2;
							continue;
						}
						nfound++;
						uint32_t e1 = adj_[v1][p1].eid; uint32_t e2 = adj_[v2][p2].eid;
						++p1; ++p2;
						//update ks
						if(k_[e1] < k-1 || k_[e2] < k-1) continue;
						kse++;
						if(k_[e1] >= k && k_[e2] >= k) {
							if(k_[e1] == k && Ins[e1] == false && ks_[e1]-- == k_[e1]) {
								q.push(e1); Ins[e1] = true;
							}
							if(k_[e2] == k && Ins[e2] == false && ks_[e2]-- == k_[e2]) {
								q.push(e2); Ins[e2] = true;
							}
						}			
					}
					else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
						++p1;
					} else {
						++p2;
					}
      			}
				ks_[e] = kse;
				k_[e]--;
				TC[k].emplace_back(e);
			}
		}//for time2triangle[t]
        for(uint32_t i = 1;i <= k_max; ++i) {
			if(TC[i].size() != kTClen[i]){
				TC_index[i].emplace_back(make_pair(t,kTClen[i]));
				kTClen[i] = TC[i].size();
			}
		}
	}
    for(uint32_t i = 0; i < m_; ++i){
		while(k_[i] > 0){
			TC[k_[i]].emplace_back(i);
			k_[i]--; 
		}  
    }

	for(uint32_t i = 1;i <= k_max; ++i) {
		if(TC[i].size() != kTClen[i]){
			TC_index[i].emplace_back(make_pair(0,kTClen[i]));
			kTClen[i] = TC[i].size();
		}
	}
	t_end = clock();
	double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC); 
    cout << "Index construction time: " << fixed << time_taken << setprecision(5); 
    cout << " sec " << endl;
	decltype(edge2times)().swap(edge2times);
	decltype(time2triangle)().swap(time2triangle);
	decltype(adj_)().swap(adj_);
	decltype(k_)().swap(k_);
	decltype(ks_)().swap(ks_);
	decltype(s_)().swap(s_);
	decltype(pos)().swap(pos);
	decltype(ord_)().swap(ord_);
	decltype(Mc)().swap(Mc);
}

void queryKDeltaTruss(uint32_t k, uint32_t t){
	uint32_t p = 0;
	while(p < TC_index[k].size() && TC_index[k][p].first > t) p++;
	if(p == TC_index[k].size()) return;
	uint32_t startP= TC_index[k][p].second;
    for(auto it = TC[k].begin()+startP; it != TC[k].end(); ++it){
    	uint32_t u = edges_[*it].first;
		uint32_t v = edges_[*it].second;
    }
}

int main(int args, char* argv[]) {
    if(args != 2) {
        cout << "! Error !" << endl;
        cout << "Usage: ./main dataset_path" << endl;
        cout << "-----------------------" << endl;
        system("ls -1 datasets");
        cout<< "-----------------------" << endl;
        exit(0);
    }
    string dataset_path = argv[1];
	tc_Index_Construction(dataset_path);
	return 0;
}
	