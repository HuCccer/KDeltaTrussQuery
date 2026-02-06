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

struct ArrayEntry{
    uint32_t vid;
    uint32_t eid;
};

// the adjacency array representation
vector<vector<ArrayEntry>> adj_;

//hashmap for invalid triangle
vector<unordered_map<uint32_t, bool>> invalidTri;

// data members
uint32_t n_;  // the # of vertices
uint32_t m_;  // the # of edges
uint32_t k_max = 0;	//the # maximum k value
// the truss numbers
VI k_;
// the edge peeling order
VI ord_;
//the position of edges in order 
VI pos;

//the timestamps of edges
vector<VI> edge2times;

time_t t_start, t_end;

//the maxId edge and minId edge of triangle
uint32_t mineId , maxeId;

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
		invalidTri.resize(m_);
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
	if (b > c){maxeId = b; mineId = a;}
	else {maxeId = c; mineId = a;}
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

void getKDeltaTruss(string dataset_path, uint32_t k, uint32_t delta){
	//Load temporal graph
	TemporalGraph tgraph = TemporalGraph(dataset_path);
	t_start = clock();
    // 1 initialization
  	// 1.1 sort the vertices in non-ascending order of degree
 	VI verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), Pred);
  	k_.resize(m_);
  	vector<vector<ArrayEntry>> A(n_);
    // 1.2 compute delta-support
	for (uint32_t v : verts) {
    	for (auto ae : adj_[v]) {
			uint32_t u = ae.vid;
			uint32_t e = ae.eid;
			if (!Pred(v, u)) continue;
			size_t pv = 0, pu = 0;
			while (pv < A[v].size() && pu < A[u].size()) {
				if (A[v][pv].vid == A[u][pu].vid) {
				uint32_t e1 = A[v][pv].eid; uint32_t e2 = A[u][pu].eid;
				uint32_t interval = get_mst(edge2times[e],edge2times[e1],edge2times[e2]); //compute minimum time span
				if(interval > delta){
                    getminmax(e, e1, e2);
                    invalidTri[mineId][maxeId] = true;//invalid triangle
                }
                else{
                    ++k_[e]; ++k_[e1]; ++k_[e2];//valid triangle
                } 
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
	// 2.decomposition
	// 2.1 sort the edges according to their supports
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
	//2.2 peeling
  	for (uint32_t i = 0; i < m_; ++i) {
        if(k_[ord_[i]] >= k) break;
    	// ASSERT(bin[k] == i);
    	const uint32_t eid = ord_[i];
    	++bin[k_[eid]];
        removeEdge(eid);
        const uint32_t v1 = edges_[eid].first;
        const uint32_t v2 = edges_[eid].second;
        size_t p1 = 0, p2 = 0;
        while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
            if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
				getminmax(eid, adj_[v1][p1].eid, adj_[v2][p2].eid);
				//the mst of triangle is greater than delta
				if(invalidTri[mineId].find(maxeId) != invalidTri[mineId].end()) {
					++p1; ++p2;
					continue;
				}
				//the mst of triangle is no greater than delta
                const uint32_t e1 = adj_[v1][p1].eid;
                const uint32_t e2 = adj_[v2][p2].eid;
                for (const uint32_t e : {e1, e2}) {
                    if (k_[e] > k_[eid]) {
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
                ++p1; ++p2;
            } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
                ++p1;
            } else {
                ++p2;
            }
        }
    }
	t_end = clock();
	cout << "-----------------------------" << endl; 
    double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC); 
    cout << "Online query time: " << fixed << time_taken << setprecision(5); 
    cout << " sec " << endl;
}

void test(){
	int edge_num = 0;
	for(int vertex = 0 ;vertex < n_ ; ++vertex){
		if(adj_[vertex].size()){
			for(auto ae : adj_[vertex]) {
				if(ae.vid > vertex) cout<<vertex<<" "<<ae.vid<<endl;
			}
			edge_num += adj_[vertex].size();
		}
	}
	cout<<"num of edges: "<<edge_num / 2<<endl;
}

void DFS(int node,vector<bool>& visited, vector<int>& component) {
    visited[node] = true;
    component.push_back(node);
    for (auto neighbor : adj_[node]) {
        if (!visited[neighbor.vid]) {
            DFS(neighbor.vid, visited, component);
        }
    }
}

void findConnectedComponents(vector<vector<int>> &components) {
    vector<bool> visited(n_, false);
    for (int i = 0; i < n_; i++) {
        if (!visited[i]) {
            vector<int> component;
            DFS(i,visited, component);
            components.push_back(component);
        }
    }
    return;
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
	uint32_t k, delta;
	cout<<"please input support constrant k (k >= 2) and time constraint delta (delta >= 0)!"<<endl;
	cin>>k>>delta;
	getKDeltaTruss(dataset_path, k - 2, delta);
}
	