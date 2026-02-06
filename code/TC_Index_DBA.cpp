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

using namespace std;

typedef vector<uint32_t> VI;
// edge type
typedef std::pair<uint32_t, uint32_t> EdgT;

// the set of edges
vector<EdgT> edges_;

// triangle type
typedef std::tuple<uint32_t, uint32_t, uint32_t> Triangle;

struct ArrayEntry{
    uint32_t vid;
    uint32_t eid;
};

// the adjacency array representation
vector<vector<ArrayEntry>> adj_;

// data members
uint32_t n_;  // the # of vertices
uint32_t m_;  // the # of edges
uint32_t k;	//the k value
uint32_t t_max; //the #maximum mts
// the support of e
vector<uint32_t> sup;


//the timestamps of edges
vector<VI> edge2times;
//delta-triangle list
vector<vector<Triangle>> time2triangle(4000);
//TC-index
vector<vector<pair<uint32_t, uint32_t>>> TC_index(100);
vector<vector<uint32_t>> TC(100);

time_t t_start, t_end;
//the maxId edge and minId edge of triangle
uint32_t mineId , maxeId;

class TemporalGraph {
    public:
    int timestamps;
    int vertices;
    int temporal_edges;

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


inline bool pred(int i, int j) {
    return adj_[i].size() > adj_[j].size() || (adj_[i].size() == adj_[j].size() && i > j);
}

uint32_t get_mst(const VI &times1,const VI &times2,const VI &times3){
	uint32_t i = 0, j = 0, k = 0;
	uint32_t n = times1.size(), m = times2.size(), p = times3.size();
	uint32_t delta_min = 20000;
	while(i < n && j < m && k < p){
		int minvalue = min({(uint32_t)times1[i], times2[j], times3[k]});
		int delta = max({(uint32_t)times1[i], times2[j], times3[k]}) - minvalue;
		if(delta < delta_min) delta_min = delta;
		if(delta_min == 0) return delta_min;
		if(minvalue == times1[i]) i++;
		else if(minvalue == times2[j]) j++;
		else k++;
	}
	return delta_min;
}

inline void getminmax(uint32_t a, uint32_t b, uint32_t c){
	if (a > b) a ^= b ^= a ^= b;
	if (a > c) a ^= c ^= a ^= c;
	if (b > c){maxeId = b; mineId = a;}
	else {maxeId = c; mineId = a;}
}

inline Triangle sortTri(uint32_t a, uint32_t b, uint32_t c){
	if (a > b) a ^= b ^= a ^= b;
	if (a > c) a ^= c ^= a ^= c;
	if (b > c) b ^= c ^= b ^= c;
	return {a, b, c};
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

//delta-decomposition
void decomph(vector<bool> &removed, VI &s_){
	queue<int> q;
	vector<int> Ins(m_,0);
	vector<unordered_map<uint32_t, bool>> invalidtri(m_);
	int cur_index = 0;
 	for(uint32_t t = t_max ; t > 0 ; --t){
		if(time2triangle[t].size() == 0) continue;
		for(auto tri:time2triangle[t]){
			uint32_t e1 = std::get<0>(tri), e2 = std::get<1>(tri), e3 = std::get<2>(tri);
			if(removed[e1] || removed[e2] || removed[e3]) continue;
			invalidtri[e1][e3] = true;
			s_[e1]--; s_[e2]--; s_[e3]--;
			
			if(s_[e1] < k && Ins[e1] == 0){
				q.push(e1);
				Ins[e1] = 1;
			}
			if(s_[e2] < k  && Ins[e2] == 0){
				q.push(e2);
				Ins[e2] = 1;
			}
			if(s_[e3] < k  && Ins[e3] == 0){
				q.push(e3);
				Ins[e3] = 1;
			}
		}
		while(!q.empty()){
			uint32_t e=q.front();
			q.pop();
			removed[e] = true;
			uint32_t supe = s_[e];
			TC[k].emplace_back(e);
			vector<std::pair<uint32_t, uint32_t>> tris; {
				const uint32_t v1 = edges_[e].first;
				const uint32_t v2 = edges_[e].second;
				size_t p1 = 0, p2 = 0;
				while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
					if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
						tris.emplace_back(adj_[v1][p1].eid, adj_[v2][p2].eid);
						++p1; ++p2;
					} else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
						++p1;
					} else {
						++p2;
					}
				}
			}
			uint32_t nfound=0;
			for (const auto tri : tris) {
				if(nfound == supe) break;
				const uint32_t e1 = tri.first;
				const uint32_t e2 = tri.second;
				if(removed[e1] || removed[e2]) continue;
				getminmax(e, e1, e2);
				if(invalidtri[mineId].find(maxeId) != invalidtri[mineId].end()) {
					continue;
				}
				else{
					invalidtri[mineId][maxeId] = true;
				}
				nfound++;
				s_[e1]--; s_[e2]--;
				if(s_[e1] < k && Ins[e1] == 0){
					q.push(e1);
					Ins[e1] = 1;
				}
				if(s_[e2] < k && Ins[e2] == 0){
					q.push(e2);
					Ins[e2] = 1;
				} 	
			}
		}
		if(cur_index != TC[k].size()){
			TC_index[k].emplace_back(make_pair(t,cur_index));
			cur_index = TC[k].size();
		}
	}
	//the case of deta = 0
	for(int i=0; i<m_; ++i){
		if(removed[i] == false){
			TC[k].emplace_back(i);
		}
	}
	if(cur_index != TC[k].size()) TC_index[k].emplace_back(make_pair(0,TC[k].size()));
}

void kDeltaDecomp(){
	// 1. initialization
	// 1.1. sort the vertices in non-ascending order of degree
	VI sup_; VI ord_; VI pos;
 	VI verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), pred);
  	// 1.2. construct delta-triangle list
  	sup_.resize(m_);
  	vector<vector<ArrayEntry>> A(n_);
	for (uint32_t v : verts) {
    	for (auto ae : adj_[v]) {
			uint32_t u = ae.vid;
			uint32_t e = ae.eid;
			if (!pred(v, u)) continue;
			size_t pv = 0, pu = 0;
			while (pv < A[v].size() && pu < A[u].size()) {
				if (A[v][pv].vid == A[u][pu].vid) {
					++sup_[A[v][pv].eid]; ++sup_[A[u][pu].eid];
					++sup_[e];
					uint32_t interval = get_mst(edge2times[e],edge2times[A[v][pv].eid],edge2times[A[u][pu].eid]);//compute minimum time span
					if(interval > t_max) t_max = interval;
					time2triangle[interval].emplace_back(sortTri(e,A[v][pv].eid,A[u][pu].eid));
					++pv; ++pu;
				} else if (pred(A[v][pv].vid, A[u][pu].vid)) {
					++pv;
				} else {
					++pu;
				}
			}
      		A[u].emplace_back(ArrayEntry{v, e});
		}
	}
	// 2. k-decomposition
	// 2.1. sort the edges according to their supports
	const uint32_t max_sup = *max_element(sup_.cbegin(), sup_.cend());
	VI bin(max_sup + 1, 0);
  	for (uint32_t eid = 0; eid < m_; ++eid) ++bin[sup_[eid]];
	for (uint32_t i = 0, start = 0; i <= max_sup; ++i) {
		start += bin[i];
		bin[i] = start - bin[i];
	}

  	ord_.resize(m_);
	pos.resize(m_);
  	for (uint32_t eid = 0; eid < m_; ++eid) {
    	pos[eid] = bin[sup_[eid]];
    	ord_[pos[eid]] = eid;
    	++bin[sup_[eid]];
  	}
	rotate(bin.rbegin(), bin.rbegin() + 1, bin.rend());
  	bin[0] = 0;
  	// 2.2. peeling
	vector<bool> removed(m_ ,false);
	vector<bool> re(m_);
	VI s_(m_);
  	for (uint32_t i = 0; i < m_; ++i) {	
    	const uint32_t eid = ord_[i];
		if(sup_[eid] != k){
			re = removed;
			s_ = sup_;
			for(int i=k+1;i<=sup_[eid];++i){
				k = i;
				//delta-decomposition
				decomph(re, s_);
			}
		}
    	++bin[sup_[eid]];
		removed[eid] = true;
    	// find triangles containing the edge with ID eid
    	vector<std::pair<uint32_t, uint32_t>> tris; {
      		const uint32_t v1 = edges_[eid].first;
      		const uint32_t v2 = edges_[eid].second;
      		size_t p1 = 0, p2 = 0;
      		while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
				if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
					tris.emplace_back(adj_[v1][p1].eid, adj_[v2][p2].eid);
					++p1; ++p2;
				} else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
					++p1;
				} else {
					++p2;
				}
      		}
    	}
		for (const auto tri : tris) {
			const uint32_t e1 = tri.first;
			const uint32_t e2 = tri.second;
			if (removed[e1] || removed[e2]) continue;
			for (const uint32_t e : {e1, e2}) {
				if (sup_[e] > k) {
				const uint32_t pe3 = bin[sup_[e]];
				const uint32_t pe = pos[e];
				if (pe3 != pe) {
					const uint32_t e3 = ord_[pe3];
					ord_[pe] = e3;
					pos[e3] = pe;
					ord_[pe3] = e;
					pos[e] = pe3;
				}
				++bin[sup_[e]];
				--sup_[e];
				}
			}
		}
  	}
}

void queryKDeltaTruss(int k, int t){
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
	TemporalGraph tgraph = TemporalGraph(dataset_path);
	t_start = clock();
	kDeltaDecomp();
	t_end = clock();
    cout << "-----------------------------" << endl; 
    double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC); 
    cout << "Index construction time: " << fixed << time_taken << setprecision(5); 
    cout << " sec " << endl;
	decltype(edge2times)().swap(edge2times);
    decltype(time2triangle)().swap(time2triangle);
    decltype(adj_)().swap(adj_);
}
	