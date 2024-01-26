#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>

#include "bits/stdc++.h"
using namespace std;
typedef long long ll;
typedef long double ld;
typedef pair<int, int> pii;
typedef vector<int> vi;

//adj list for connected components
typedef vector<vector<int>> ADJLIST;

#define FOR(i, hi) for(int i=0; i<hi; ++i)
#define FORR(i, hi) for(int i=hi-1; i>=0; --i)
#define For(i, lo, hi) for(int i=lo; i<hi; ++i)
#define Forr(i, hi, lo) for(int i=hi-1; i>=lo; --i)

struct PepMetadata {
    double pepmass;
    int global_scan;
};

struct Match {
    int a, b;
    double massprod;
    bool shared;
};

struct SpectrumInfo {
    int ind;
    int scan;
    double true_mz;
    double mass;
};

struct PeakInfo {
    int mz_bin;
    double true_mz;
    double mass;
};

struct ProductInfo {
    int refscan;
    double product;
    double product_shared;
    double product_shifted;
};

chrono::duration<double> tot_sort;
typedef vector<vector<SpectrumInfo>> mztospectra_t;
//typedef unordered_map<int, vector<PeakInfo>> spectratomz_t;
typedef vector<vector<PeakInfo>> spectratomz_t;
typedef pair<mztospectra_t, spectratomz_t> table_t;
static double TOLERANCE = 0.01;
static const int TOPK = 10;
static const int MIN_MATCHED_PEAKS = 0;
static const double BAD_THRESH = 0.1;
static double THRESHOLD = 0.7;
static int NUM_FILES_MAX = 50000;
static double SELECTIONRANGE = 200.0;
static int FILTERK = 20;
static int TOPProducts = 10;
static int BruteForceThresh = 1000;
static int MINMATCHES = 6;
static int MINPEAKS = 6;
int shifted_offset = 6000;


//code for calculating connected components using BFS
vector<int> BFS_helper(ADJLIST &adj_list, int v, vector<bool> &discovered, vector<int> &components_info, int component_index) {
    // create a queue for doing BFS
    vector<int> result;
    queue<int> q;
    // mark the source vertex as discovered
    discovered[v] = true;
 
    // enqueue source vertex
    q.push(v);
    components_info[v] = component_index; 
    result.push_back(v);
 
    // loop till queue is empty
    while (!q.empty()) {
        // dequeue front node and print it
        v = q.front();
        q.pop();

        // do for every edge (v, u)
        for (int u: adj_list[v]) {
            if (!discovered[u]) {
                // mark it as discovered and enqueue it
                discovered[u] = true;
                q.push(u);
                components_info[u] = component_index; 
                result.push_back(u);
            }
        }
    }
    return result;
}


int calculate_degree(ADJLIST &adj_list, vector<int> &component_vertices) {
    int res = 0;
    for (int v: component_vertices) {
        res += adj_list[v].size();
    }
    return res / 2;
}

vector<int> BFS_all(ADJLIST &adj_list, int max_node, string connected_components_out) {
    ofstream cluster_info_out_write(connected_components_out);
    cluster_info_out_write << "connected component index " << "\t" << "component size" << "\t" <<" average_degree " << endl;
    vector<bool> visited;
    vector<int> components_info;
    for (int v = 0; v < max_node; v++) {
        visited.push_back(false);
        components_info.push_back(-1);
    }
    int component_index = 0;
    for (int node = 0; node < max_node; node ++) {
        if (visited[node]) {
            continue;
        } else {
            vector<int> component_vertex_list = BFS_helper(adj_list, node, visited, components_info, component_index);
            int component_deg = calculate_degree(adj_list, component_vertex_list);
            int component_size = component_vertex_list.size();
            cluster_info_out_write << component_index << "\t" << component_size << "\t" << component_deg << endl;
            component_index += 1;
        }
    }
    cluster_info_out_write.close();
    return components_info;
}



//code for calculating connected components
void DFS_helper(int node, int component_index, ADJLIST &adj_list, vector<bool> &visited, vector<int> &components_info) {
    if (visited[node]) {
        return;
    } else {
        visited[node] = true;
        components_info[node] = component_index;
        for (int j = 0; j< adj_list[node].size(); j++) {
            int neighbor = adj_list[node][j];
            DFS_helper(neighbor, component_index, adj_list, visited, components_info);
        }
    }
}


vector<int> DFS_all(ADJLIST &adj_list, int max_node) {
    vector<bool> visited;
    vector<int> components_info;
    for (int v = 0; v < max_node; v++) {
        visited.push_back(false);
        components_info.push_back(-1);
    }
    int component_index = 0;
    for (int node = 0; node < max_node; node ++) {
        if (visited[node]) {
            continue;
        } else {
            DFS_helper(node, component_index, adj_list, visited, components_info);
            component_index += 1;
        }
    }
    return components_info;
}





vector<string> get_files(string filename) {
    vector<string> filenames;

    ifstream f(filename);
    string a;

    while(f >> a && filenames.size() < NUM_FILES_MAX) {
		cout << a << endl;
        filenames.push_back(a);
    }
    f.close();

    return filenames;
}

bool sort_tuple(const tuple<int, double, double>& a, const tuple<int, double, double>& b) {
    return (get<2>(a) > get<2>(b));
}

bool examine_filenames(vector<string> &filenames) {
    bool all_file_good = true;
    for (string filename : filenames) {
        ifstream infile(filename);
        if (!infile.good()) {
            cout << "file not exist: " << filename << endl;
            all_file_good = false;
        }
        infile.close();
    }
    return all_file_good;
}


pair<pair<table_t, table_t>, vector<PepMetadata>> parse_inputs_from_files(vector<string> &filenames) {
    pair<pair<table_t, table_t>, vector<PepMetadata>> input;
    table_t& table_shared = input.first.first;
    table_t& table_shifted = input.first.second;

    mztospectra_t& mztospectra_shared = table_shared.first;
    mztospectra_t& mztospectra_shifted = table_shifted.first;

    spectratomz_t& spectratomz_shared = table_shared.second;
    spectratomz_t& spectratomz_shifted = table_shifted.second;

    vector<PepMetadata>& metadata = input.second;

    mztospectra_shared.resize(200000);
    mztospectra_shifted.resize(200000);

    double increase_val = 0.0;
    unsigned int total_peaks = 0;
    int spectra_counter = 0;
    int empty_counter = 0;
    int scan_tracer = 0;
    int global_scan_tracer = 0;
    for(string filename : filenames) {
		cout << filename << endl;
        ifstream f(filename);
        string a;
        int max_mz_shared = 0;
        int max_mz_shifted = 0;
        while(getline(f, a)) {
            if(a == "BEGIN IONS") {
                double pepmass = -1;
                //int scan = -1;

                while(getline(f, a)) {
                    if(a.rfind("PEPMASS", 0) == 0)
                        pepmass = stod(a.substr(8));

                    if(a.rfind("RTINSECONDS",0) == 0) {
                        //scan = stoi(a.substr(6));

                        if(pepmass == -1) break;

                        vector<vector<tuple<int, double, double>>> curv_shared;
                        vector<vector<tuple<int, double, double>>> curv_shifted;
                        long double magnitude = 0;
                        bool nonempty = false;
                        while(getline(f, a)) {
                            if(a == "END IONS") break;

                            double mz, mass;
                            stringstream ss(a);
                            ss >> mz >> mass;

                            total_peaks += 1;
                            int ind_shared = round(mz / TOLERANCE);
                            int ind_shifted = round((pepmass - mz + shifted_offset) / TOLERANCE);
                            
                            if(ind_shifted < 0)
                                increase_val = max(increase_val, -(pepmass - mz + shifted_offset));

                            int filter_k_index_shared = round(max_mz_shared / SELECTIONRANGE);

                            if (filter_k_index_shared >= curv_shared.size()) {
                                curv_shared.resize(filter_k_index_shared * 3 / 2 + 1);
                            }
                            if (filter_k_index_shared >= curv_shifted.size()) {
                                curv_shifted.resize(filter_k_index_shared * 3 / 2 + 1);
                            }
                            tuple<int, double, double> elem1{ind_shared, mz, mass};
                            curv_shared[filter_k_index_shared].push_back(elem1);
                            tuple<int, double, double> elem2{ind_shifted, pepmass - mz + shifted_offset, mass};
                            curv_shifted[filter_k_index_shared].push_back(elem2);

                            nonempty = true;
                        }
                        vector<tuple<int, double, double>> filtered_shared_peaks;
                        vector<tuple<int, double, double>> filtered_shifted_peaks;

                        for (int filter_index = 0; filter_index < curv_shared.size(); filter_index ++) {
                            vector<tuple<int, double, double>> primitive_bin_shared = curv_shared[filter_index];
                            vector<tuple<int, double, double>> primitive_bin_shifted = curv_shifted[filter_index];
                            
                            if (primitive_bin_shared.size() >= FILTERK) {
                                sort(primitive_bin_shared.begin(), primitive_bin_shared.end(), sort_tuple);
                                sort(primitive_bin_shifted.begin(), primitive_bin_shifted.end(), sort_tuple);
                                curv_shared[filter_index] = {primitive_bin_shared.begin(), primitive_bin_shared.begin()+FILTERK};
                                curv_shifted[filter_index] = {primitive_bin_shifted.begin(), primitive_bin_shifted.begin()+FILTERK};
                                for (int j = 0; j < FILTERK; j++) {
                                    filtered_shared_peaks.push_back(primitive_bin_shared[j]);
                                    filtered_shifted_peaks.push_back(primitive_bin_shifted[j]);
                                    double mass = get<2>(primitive_bin_shared[j]);
                                    magnitude += mass * mass;
                                }
                            } else {
                                for (int j = 0; j < primitive_bin_shared.size(); j++) {
                                    filtered_shared_peaks.push_back(primitive_bin_shared[j]);
                                    filtered_shifted_peaks.push_back(primitive_bin_shifted[j]);
                                    double mass = get<2>(primitive_bin_shared[j]);
                                    magnitude += mass * mass;
                                } 
                            }
                        }

                        if((nonempty && increase_val == 0.0) && (magnitude > 0)) {
                            magnitude = sqrt(magnitude);
                            int scan = scan_tracer;
                            scan_tracer += 1;

							if (scan >= spectratomz_shared.size()) {
								spectratomz_shared.resize(scan + 1);
								spectratomz_shifted.resize(scan + 1);
								metadata.resize(scan + 1);
							}
                            
                            for(int i = 0; i < filtered_shared_peaks.size(); ++i) {
                                tuple<int, double, double> p = filtered_shared_peaks[i];
                                int ind_shared = get<0>(p);
                                double mz = get<1>(p);
                                double mass = get<2>(p);
								max_mz_shared = max(max_mz_shared, ind_shared);
                                spectratomz_shared[scan].push_back({ind_shared, mz, (double)(mass / magnitude)});
                                if(max_mz_shared >= mztospectra_shared.size()) {
                                    mztospectra_shared.resize(max_mz_shared * 3 / 2 + 1);
                                }
                                mztospectra_shared[ind_shared].push_back({i, scan, mz, (double)(mass / magnitude)});
                            }
                            
                            for(int i = 0; i < filtered_shifted_peaks.size(); ++i) {
                                tuple<int, double, double> p = filtered_shifted_peaks[i];
                                int ind_shifted = get<0>(p);
                                double mz_shifted = get<1>(p);
                                double mass = get<2>(p);
                                max_mz_shifted = max(max_mz_shifted, ind_shifted);
                                spectratomz_shifted[scan].push_back({ind_shifted, mz_shifted, (double)(mass / magnitude)});
                                if(max_mz_shifted >= mztospectra_shifted.size()) {
                                    mztospectra_shifted.resize(max_mz_shifted * 3 / 2 + 1);
                                }
                                mztospectra_shifted[ind_shifted].push_back({i, scan, mz_shifted, (double)(mass / magnitude)});
                            }
                            
                            metadata[scan] = PepMetadata{pepmass, global_scan_tracer};
                            global_scan_tracer += 1;
                            spectra_counter += 1;

                        } else {
                            global_scan_tracer += 1;
                            empty_counter += 1;
                        }

                        break;
                    }
                    if(a == "END IONS") break;
                }
            }
        }

        f.close();
    }
    
    if(increase_val > 0.0) {
        cout << "Increase shifted offset by at least " << increase_val << endl;
        cout << "minimum shifted offset required " << increase_val + shifted_offset << endl;
        assert(false);
    }
    cout << "total number of peaks: " << total_peaks << endl;
    cout << "total number of non-empty spectras in data source: " << spectra_counter << endl;
    cout << "total number of empty spectras in data source: " << empty_counter << endl;
    return input;
}


struct scan_comp {
    inline bool operator() (const SpectrumInfo& s1, const SpectrumInfo& s2) {
        return (s1.scan < s2.scan);
    }
    inline bool operator() (const SpectrumInfo& s1, double scan) {
        return (s1.scan < scan);
    }
    inline bool operator() (double scan, const SpectrumInfo& s2) {
        return (scan < s2.scan);
    }
};

/*
void update_matches(mztospectra_t& mztospectra, int queryscan, int curi, int bin, double mz, double mass, unordered_map<int, vector<Match>>& best_matches) {
    
    auto t1 = upper_bound(mztospectra[bin].begin(), mztospectra[bin].end(), queryscan, scan_comp());
    for(; t1 != mztospectra[bin].end(); ++t1) {
        best_matches[t1->scan].push_back({curi, t1->ind, mass * t1->mass});
    }
    
    if(bin > 0) {
        auto t2 = upper_bound(mztospectra[bin-1].begin(), mztospectra[bin-1].end(), queryscan, scan_comp());
        for(; t2 != mztospectra[bin-1].end(); ++t2) {
            if(abs(mz - t2->true_mz) <= TOLERANCE)
                best_matches[t2->scan].push_back({curi, t2->ind, mass * t2->mass});
        }
    }
    if(bin < mztospectra.size()-1) {
        auto t2 = upper_bound(mztospectra[bin+1].begin(), mztospectra[bin+1].end(), queryscan, scan_comp());
        for(; t2 != mztospectra[bin+1].end(); ++t2) {
            if(abs(mz - t2->true_mz) <= TOLERANCE)
                best_matches[t2->scan].push_back({curi, t2->ind, mass * t2->mass});
        }
    }
}
*/


struct pair_comp {
    inline bool operator() (const Match& p1, const Match& p2) {
        return (p1.massprod > p2.massprod);
    }
};

bool cmp(Match x, Match y) {
	return (x.massprod > y.massprod);
}


vector<ProductInfo> brute_force_selection_storage(vector<ProductInfo> &product_storage, int max_selection) {
    vector<ProductInfo> result;
    set<int> seen;
    int ranking = 0;
    int product_storage_size = product_storage.size();
    int total_preserved = min(max_selection, product_storage_size);
    assert(total_preserved <= product_storage.size());
    for (int iters = 0; iters < total_preserved; iters++) {
        double max_val = 0.0;
        int max_idx = -1;
        for (int i = 0; i < product_storage.size(); i++) {
            if (seen.find(i) == seen.end()) {
                double prod = product_storage[i].product;
                assert(prod >= THRESHOLD);
                if (prod > max_val) {
                    max_val = prod;
                    max_idx = i;
                }
            }
        }

        //cout << "product storage size: " << product_storage.size() << ", and there are " << seen.size() << " values "<< endl;
        assert(max_idx >= 0);
        assert(max_idx < product_storage.size());
        seen.insert(max_idx);
        result.push_back({product_storage[max_idx].refscan, product_storage[max_idx].product, product_storage[max_idx].product_shared, product_storage[max_idx].product_shifted});
    }
    return result;
}

bool cmp_prods(ProductInfo x, ProductInfo y) {
	return (x.product > y.product);
}

vector<ProductInfo> quick_sort_selection_storage(vector<ProductInfo> &product_storage, int max_selection) {
    sort(product_storage.begin(), product_storage.end(), cmp_prods);
    int product_storage_size = product_storage.size();
    int total_preserved = min(max_selection, product_storage_size);
    vector<ProductInfo> result = {product_storage.begin(), product_storage.begin() + total_preserved};
    return result;
}





unordered_map<int, unordered_map<int, double>> all_versus_all(pair<table_t, table_t>& input, vector<PepMetadata>& masses, string output_file, string cluster_info_in, string connected_components_out) {
    unordered_map<int, unordered_map<int, double>> matches;

    table_t& table_shared = input.first;
    table_t& table_shifted = input.second;
    mztospectra_t& mztospectra_shared = table_shared.first;
    mztospectra_t& mztospectra_shifted = table_shifted.first;

    spectratomz_t& spectratomz_shared = table_shared.second;
    spectratomz_t& spectratomz_shifted = table_shifted.second;
    
    int total_specnum = spectratomz_shared.size();
    int spec_counter = 0;
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    int scan_tracer = 0;
    unordered_map<int, vector<ProductInfo>> all_results;
    for(int scan = 0; scan < spectratomz_shared.size(); scan ++) {
        vector<double> upperbound_prod;
        vector<int> match_counter;

        upperbound_prod.resize(spectratomz_shared.size());
        match_counter.resize(spectratomz_shared.size());
        
        vector<ProductInfo> product_storage;
                

        assert(product_storage.size() == 0);
        unordered_map<int, vector<Match>> best_matches;
        vector<PeakInfo>& peaks_share = spectratomz_shared[scan];
        vector<PeakInfo>& peaks_shift = spectratomz_shifted[scan];
        for(int i = 0; i < peaks_share.size(); ++i) {
            PeakInfo peak_shared = peaks_share[i];
            PeakInfo peak_shifted = peaks_shift[i];

            int bin_shared = peak_shared.mz_bin;
            int bin_shifted = peak_shifted.mz_bin;

            double mz_shared = peak_shared.true_mz;
            double mz_shifted = peak_shifted.true_mz;

            double mass = peak_shared.mass;
            if (peak_shifted.mass != mass) {
                cout << "peak shifted: " << peak_shifted.mass << endl;
                cout << "peak shared: " << peak_shared.mass << endl;
            }
            assert(peak_shifted.mass == mass);


            auto t1 = upper_bound(mztospectra_shared[bin_shared].begin(), mztospectra_shared[bin_shared].end(), scan, scan_comp());
            for(; t1 != mztospectra_shared[bin_shared].end(); ++t1) {
                upperbound_prod[t1->scan] += mass * t1->mass;
                match_counter[t1->scan] += 1;
                
            }    
            if(bin_shared > 0) {
                auto t2 = upper_bound(mztospectra_shared[bin_shared-1].begin(), mztospectra_shared[bin_shared-1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shared[bin_shared-1].end(); ++t2) {
                    if(abs(mz_shared - t2->true_mz) <= TOLERANCE) {
                        upperbound_prod[t2->scan] += mass * t2->mass;
                        match_counter[t2->scan] += 1;
                    }
                }
            }
            if(bin_shared < mztospectra_shared.size()-1) {
                auto t2 = upper_bound(mztospectra_shared[bin_shared+1].begin(), mztospectra_shared[bin_shared+1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shared[bin_shared+1].end(); ++t2) {
                    if(abs(mz_shared - t2->true_mz) <= TOLERANCE) {
                        upperbound_prod[t2->scan] += mass * t2->mass;
                        match_counter[t2->scan] += 1;
                    }
                }
            }
            auto t3 = upper_bound(mztospectra_shifted[bin_shifted].begin(), mztospectra_shifted[bin_shifted].end(), scan, scan_comp());
            for(; t3 != mztospectra_shifted[bin_shifted].end(); ++t3) {
                upperbound_prod[t3->scan] += mass * t3->mass;
                match_counter[t3->scan] += 1;
            }
                
            if(bin_shifted > 0) {
                auto t2 = upper_bound(mztospectra_shifted[bin_shifted-1].begin(), mztospectra_shifted[bin_shifted-1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shifted[bin_shifted-1].end(); ++t2) {
                    if(abs(mz_shifted - t2->true_mz) <= TOLERANCE) {
                        upperbound_prod[t2->scan] += mass * t2->mass;
                        match_counter[t2->scan] += 1;
                    }
                }
            }
            if(bin_shifted < mztospectra_shifted.size()-1) {
                auto t2 = upper_bound(mztospectra_shifted[bin_shifted+1].begin(), mztospectra_shifted[bin_shifted+1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shifted[bin_shifted+1].end(); ++t2) {
                    if(abs(mz_shifted - t2->true_mz) <= TOLERANCE) {
                        upperbound_prod[t2->scan] += mass * t2->mass;
                        match_counter[t2 -> scan] += 1;
                    }
                }
            }
        }

        for(int i = 0; i < peaks_share.size(); ++i) {
            PeakInfo peak_shared = peaks_share[i];
            PeakInfo peak_shifted = peaks_shift[i];

            int bin_shared = peak_shared.mz_bin;
            int bin_shifted = peak_shifted.mz_bin;

            double mz_shared = peak_shared.true_mz;
            double mz_shifted = peak_shifted.true_mz;

            double mass = peak_shared.mass;
            assert(peak_shifted.mass == mass);


            auto t1 = upper_bound(mztospectra_shared[bin_shared].begin(), mztospectra_shared[bin_shared].end(), scan, scan_comp());
            for(; t1 != mztospectra_shared[bin_shared].end(); ++t1) {
                if ((upperbound_prod[t1 -> scan] < THRESHOLD) || match_counter[t1->scan] < MINMATCHES) {
                    continue;
                }
                best_matches[t1->scan].push_back({i, t1->ind, mass * t1->mass, true});
            }    
            if(bin_shared > 0) {
                auto t2 = upper_bound(mztospectra_shared[bin_shared-1].begin(), mztospectra_shared[bin_shared-1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shared[bin_shared-1].end(); ++t2) {
                    if ((upperbound_prod[t2 -> scan] < THRESHOLD) || match_counter[t2->scan] < MINMATCHES) {
                        continue;
                    }
                    if(abs(mz_shared - t2->true_mz) <= TOLERANCE) {
                        best_matches[t2->scan].push_back({i, t2->ind, mass * t2->mass, true});
                    }
                }
            }
            if(bin_shared < mztospectra_shared.size()-1) {
                auto t2 = upper_bound(mztospectra_shared[bin_shared+1].begin(), mztospectra_shared[bin_shared+1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shared[bin_shared+1].end(); ++t2) {
                    if ((upperbound_prod[t2 -> scan] < THRESHOLD) || match_counter[t2->scan] < MINMATCHES) {
                        continue;
                    }
                    if(abs(mz_shared - t2->true_mz) <= TOLERANCE) {
                        best_matches[t2->scan].push_back({i, t2->ind, mass * t2->mass, true});
                    }
                }
            }
            auto t3 = upper_bound(mztospectra_shifted[bin_shifted].begin(), mztospectra_shifted[bin_shifted].end(), scan, scan_comp());
            for(; t3 != mztospectra_shifted[bin_shifted].end(); ++t3) {
                if ((upperbound_prod[t3 -> scan] < THRESHOLD) || match_counter[t3->scan] < MINMATCHES) {
                    continue;
                }
                best_matches[t3->scan].push_back({i, t3->ind, mass * t3->mass, false});
            }
                
            if(bin_shifted > 0) {
                auto t2 = upper_bound(mztospectra_shifted[bin_shifted-1].begin(), mztospectra_shifted[bin_shifted-1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shifted[bin_shifted-1].end(); ++t2) {
                    if ((upperbound_prod[t2 -> scan] < THRESHOLD) || match_counter[t2->scan] < MINMATCHES) {
                        continue;
                    }
                    if(abs(mz_shifted - t2->true_mz) <= TOLERANCE) {
                        best_matches[t2->scan].push_back({i, t2->ind, mass * t2->mass, false});
                    }
                }
            }
            if(bin_shifted < mztospectra_shifted.size()-1) {
                auto t2 = upper_bound(mztospectra_shifted[bin_shifted+1].begin(), mztospectra_shifted[bin_shifted+1].end(), scan, scan_comp());
                for(; t2 != mztospectra_shifted[bin_shifted+1].end(); ++t2) {
                    if ((upperbound_prod[t2 -> scan] < THRESHOLD) || match_counter[t2->scan] < MINMATCHES) {
                        continue;
                    }
                    if(abs(mz_shifted - t2->true_mz) <= TOLERANCE) {
                        best_matches[t2->scan].push_back({i, t2->ind, mass * t2->mass, false});
                    }
                }
            }
        }
        for(auto it = best_matches.begin(); it != best_matches.end(); ++it) {
			set<int> seen_a, seen_b;
            int refscan = it->first;
            if (refscan >= masses.size()) {
                cout << "invalid refscan appeared when running scan number " << scan << " invalid refscan number is " << refscan << endl;
                cout << masses[refscan].pepmass << endl;
            }
            assert(refscan < masses.size());	
            if (upperbound_prod[refscan] < THRESHOLD) {
                continue;
            }
			
            vector<Match>& prods = it->second;

            int querypeaks = (int) peaks_share.size();
			int refpeaks = (int) spectratomz_shared[refscan].size();
            int max_matches = min(TOPK, min(querypeaks, refpeaks));

            int num_found = 0;
			
            double tmp_res = 0.0;
            double shared_part = 0.0;
            double shifted_part = 0.0;
            sort(prods.begin(), prods.end(), cmp);
			for (auto match_iter = prods.begin(); match_iter!= prods.end(); ++match_iter) {
			    if (num_found >=  max_matches) {
					break;
				}
				Match m = *match_iter;
				if (seen_a.find(m.a) == seen_a.end() && seen_b.find(m.b) == seen_b.end()) {
					tmp_res += m.massprod;
					seen_a.insert(m.a);
					seen_b.insert(m.b);
                    if (m.shared == true) {
                        shared_part += m.massprod;
                    } else {
                        shifted_part += m.massprod;
                    }
					num_found ++;
                }
            }
            if (tmp_res >= THRESHOLD) {
                product_storage.push_back({refscan, tmp_res, shared_part, shifted_part});
            }
        }
        
        if (product_storage.size() > BruteForceThresh) {
            all_results[scan] = brute_force_selection_storage(product_storage, TOPProducts);
        } else {
            all_results[scan] = quick_sort_selection_storage(product_storage, TOPProducts);
        }
        
        //product_storage = quick_sort_selection_storage(product_storage, TOPProducts);
        /*
        for(auto it = M.begin(); it != M.end();) {
            if(it->second < THRESHOLD) M.erase(it++);
            else ++it;
        }
        */
        //all_results[scan] = product_storage;
        spec_counter += 1;
        if (spec_counter % 100000 == 0) {
          end = chrono::steady_clock::now();
          chrono::duration<double> elapsed_seconds = end-start; 
          cout << "processed " << spec_counter << " number of spectras; Taking " << elapsed_seconds.count() << " seconds so far" << endl;
        }
    }
    /*
    cout << "begin postprocessing" << endl;
    auto writing_start = chrono::steady_clock::now();
    //calculating a connected component graph

    cout << "generating connected components" << endl;
    vector<vector<int>> adj_list;
    adj_list.resize(spectratomz_shared.size());
    for(auto iter1 = all_results.begin(); iter1 != all_results.end(); iter1++) {
        int scan = iter1 -> first;
        vector<ProductInfo> &product_storage = iter1 -> second;
        for (int k = 0; k < product_storage.size(); k++) {
            int refscan = product_storage[k].refscan;
            adj_list[scan].push_back(refscan);
            adj_list[refscan].push_back(scan); 
        }
    }
    cout << "starting to run DFS" << endl;
    vector<int> connected_components = DFS_all(adj_list, spectratomz_shared.size());

    */

    cout << "start writing ouput" << endl;
    ofstream intermediate_prints;
    intermediate_prints.open(output_file);
    int prods_counter = 0;
    //intermediate_prints << "scan_1" << "\t" << "mz_1" << "\t" << "scan_2" << "\t" << "mz_2" << "\t" << "dot_product" << "\t" << "dot_product_shared" << "\t" << "dot_product_shifted"<< "\t" << "connected component index"<< endl;
    intermediate_prints << "scan_1" << "\t" << "mz_1" << "\t" << "scan_2" << "\t" << "mz_2" << "\t" << "dot_product" << "\t" << "dot_product_shared" << "\t" << "dot_product_shifted"<< endl;
    for(auto iter1 = all_results.begin(); iter1 != all_results.end(); iter1++) {
        int scan = iter1 -> first;
        //cout <<"begin writing scan :"<<scan << endl;
        vector<ProductInfo> &product_storage = iter1 -> second;
        double pepmass1 = masses[scan].pepmass;
        for (int k = 0; k < product_storage.size(); k++) {
            int refscan = product_storage[k].refscan;
            double product = product_storage[k].product;
            double product_shared = product_storage[k].product_shared;
            double product_shifted = product_storage[k].product_shifted;
            double pepmass2 = masses[refscan].pepmass;

            double global_refscan = masses[refscan].global_scan;
            double global_scan = masses[scan].global_scan;
            //assert(connected_components[scan] == connected_components[refscan]);
            //int component = connected_components[scan];
            intermediate_prints << global_scan << "\t" << pepmass1 << "\t" << global_refscan << "\t" << pepmass2 << "\t" << product << "\t" << product_shared << "\t" << product_shifted << endl;
            prods_counter ++;
        }
    }
    intermediate_prints.close();

    cout << "begin postprocessing" << endl;
    auto writing_start = chrono::steady_clock::now();
    /*calculating a connected component graph*/

    cout << "generating connected components" << endl;
    vector<vector<int>> adj_list;
    adj_list.resize(spectratomz_shared.size());
    for(auto iter1 = all_results.begin(); iter1 != all_results.end(); iter1++) {
        int scan = iter1 -> first;
        vector<ProductInfo> &product_storage = iter1 -> second;
        for (int k = 0; k < product_storage.size(); k++) {
            int refscan = product_storage[k].refscan;
            adj_list[scan].push_back(refscan);
            adj_list[refscan].push_back(scan); 
        }
    }

    vector<int> component_capacity;
    component_capacity.resize(10000);
    int max_component_index = 0;
    vector<int> connected_components = BFS_all(adj_list, spectratomz_shared.size(), connected_components_out);
    for (int scan = 0; scan < connected_components.size(); scan++) {
        int component_index = connected_components[scan];
        if (component_index > max_component_index) {
            max_component_index = component_index;
        }
        if (component_index >= component_capacity.size()) {
            component_capacity.resize(component_index * 2 + 1);
        }
        component_capacity[component_index] += 1;
    }

    sort(component_capacity.begin(), component_capacity.end());
    
    cout << " component sizes " << "\t" << " component counter " << endl;
    int size_tracer = 0;
    int size_counter = 0;
    for (int k = 0; k < component_capacity.size(); k++) {
        if (component_capacity[k] == size_tracer) {
            size_counter += 1;
        } else {
            if (size_tracer > 0) {
                cout << size_tracer << "\t" << size_counter << endl;
            }
            while (size_tracer < component_capacity[k]) {
                size_tracer += 1;
            }
            assert(size_tracer == component_capacity[k]);
            size_counter = 1;
        }
    }



    /*
    cout << "starting to run DFS" << endl;
    vector<int> connected_components = DFS_all(adj_list, spectratomz_shared.size());
    vector<int> component_capacity;
    component_capacity.resize(10000);
    int max_component_index = 0;
    ofstream cluster_info_out_write(connected_components_out);
    for (int scan = 0; scan < connected_components.size(); scan++) {
        int component_index = connected_components[scan];
        if (component_index > max_component_index) {
            max_component_index = component_index;
        }
        if (component_index >= component_capacity.size()) {
            component_capacity.resize(component_index * 2);
        }
        component_capacity[component_index] += 1;
    }
    cluster_info_out_write << "connected component index " << "\t" << "component size" << endl;
    for (int component_idx = 0; component_idx < max_component_index + 1; component_idx += 1) {
        cluster_info_out_write << component_idx << "\t" << component_capacity[component_idx] << endl;
    }
    cluster_info_out_write.close(); 
    */

    /*
    ifstream cluster_info_in_read(cluster_info_in);
    ofstream cluster_info_out_write(connected_components_out);
    string cluster_line;
    int line_counter = 0;
    int query_scan = 0;
    while(getline(cluster_info_in_read, cluster_line)) {
        if (line_counter == 0) {
            cluster_info_out_write << "cluster_idx" << "\t" << "average mz"<< "\t" << "average RT" << "\t" << "intensity sum" << "\t" << "num spectra" << "\t" << "connected component index"<< endl;
            line_counter += 1;
            continue;
        }
        stringstream ss(cluster_line);
        int cluster_idx, num_spectra;
        double pepmass_sum_avg, RT_avg, intensity_sum; 
        ss >> cluster_idx >> pepmass_sum_avg >> RT_avg >> intensity_sum >> num_spectra;
        if (query_scan >= masses.size()) {
            break;
        } else {
            int global_scan = line_counter - 1;
            int query_global_scan = masses[query_scan].global_scan;
            if (global_scan < query_global_scan) {
                continue;
            } else {
                assert(global_scan == query_global_scan);
                int component_idx = connected_components[query_scan];
                cluster_info_out_write << cluster_idx << "\t" << pepmass_sum_avg << "\t" << RT_avg << "\t" << intensity_sum << "\t" << num_spectra << "\t" << component_idx << endl;
                query_scan += 1;
            }
        }
        line_counter += 1;
    }
    cluster_info_in_read.close();
    cluster_info_out_write.close();
    */

    auto writing_end = chrono::steady_clock::now();
    chrono::duration<double> writing_seconds = end-start; 
    cout <<  "writing " << prods_counter << " number of products; Taking " << writing_seconds.count() << " seconds so far " << endl;
    return matches;
}

int main(int argc, char* argv[]) {
    string input_type;
    string output_type;
    string input_file;
    string output_file;
    string threshold_catcher;
    string cluster_info_in = "cluster_info.tsv";
    string connected_components_out = "cluster_info_with_cc.tsv";
    string components_output_file = "components.tsv";
    string arguments = "";
    for (int k = 0; k < argc; k++) {
        arguments.append(argv[k]);
        arguments.append(" ");
    }

    cout << arguments << endl;
    bool single_file = true;
    if (argc < 5) {
        cout << "instructions for running the code: " << endl;
        cout << "usage: ./score [-i/-l <input_path>] -o <output_path> -t <default = 0.7> --peak_resolution <default=0.01>" << endl;
        cout << "\n"<< endl;
        cout << "--help: print instructions"<< endl;
        cout << "-i <input_path> : calculate pairwise dot product between spectra from a single input file " << endl; 
        cout << "-l <input_path> : calculate pairwise dot product between spectra from a list of input files" << endl; 
        cout << "-o <output_path>: output path for dot product results" << endl;
        cout << "--peak_resolution <default=0.01> (optional): the minimum mz difference for adding the product between two peaks in calculating the spectra pairwise dotproduct"<< endl;
        cout << "-t <default=0.7> (optional): the minimum threshold of preserving dot product between two spectra"<< endl;
        cout << "-f X<default=200.0> Y<default=20> (optional): add filtering to peaks in the input spectra. X: the range of filtering peaks; Y: preserving the top Y peaks in each mz range X during filtering" << endl;
        cout << "--min_peakmatches <default=6> (optional): the minimum number of matches between two spectras' shifted and shared peak locations to keep their dot product" << endl;
        cout << "--shifted_offset <default=6000> (optional): the shifted offset for calculating shifted peak mz values " << endl;
        cout << "--edge_limit <default = 10> (optional): preserving only top N products between query spectra and refscans" << endl;
        cout << "--cluster_info_in <default = cluster_info.tsv> (optional) information of the clustered results" << endl;
        cout << "--connected_components_out <default = cluster_info_with_cc.tsv> (optional) clustered results with connected components information" << endl; 
        cout << "\n" << endl;
        return 0; 
    }
    for (int i = 1; i < argc; i++) {
        //cout << argc << " " << argv[i] << endl;
        if ((i == argc - 1) && (i==1)) {
            //cout << argc << " " << argv[i] << endl;
            if (strcmp(argv[i], "--help") == 0) {
                cout << "instructions for running the code: " << endl;
                cout << "usage: ./score [-i/-l <input_path>] -o <output_path> -t <default = 0.7> --peak_resolution <default=0.01>" << endl;
                cout << "\n"<< endl;
                cout << "--help: print instructions"<< endl;
                cout << "-i <input_path> : calculate pairwise dot product between spectra from a single input file " << endl; 
                cout << "-l <input_path> : calculate pairwise dot product between spectra from a list of input files" << endl; 
                cout << "-o <output_path>: output path for dot product results" << endl;
                cout << "--peak_resolution <default=0.01> (optional): the minimum mz difference for adding the product between two peaks in calculating the spectra pairwise dotproduct"<< endl;
                cout << "-t <default=0.7> (optional): the minimum threshold of preserving dot product between two spectra"<< endl;
                cout << "-f X<default=200.0> Y<default=20> (optional): add filtering to peaks in the input spectra. X: the range of filtering peaks; Y: preserving the top Y peaks in each mz range X during filtering" << endl;
                cout << "--min_peakmatches <default=6> (optional): the minimum number of matches between two spectras' shifted and shared peak locations to keep their dot product" << endl;
                cout << "--shifted_offset <default=6000> (optional): the shifted offset for calculating shifted peak mz values " << endl;
                cout << "--edge_limit <default = 10> (optional): preserving only top N products between query spectra and refscans" << endl;
                cout << "--cluster_info_in <default = cluster_info.tsv> (optional) information of the clustered results" << endl;
                cout << "--connected_components_out <default = cluster_info_with_cc.tsv> (optional) clustered results with connected components information" << endl; 
                cout << "\n" << endl;
                return 0; 
            }
        }
        if (i + 1 < argc) {
            if ((strcmp(argv[i], "-i" ) == 0) && (input_type == "")) {                 
                input_type = "-i";    // The next value in the array is your value
                input_file = argv[i+1];
                single_file = true;
                i++;    // Move to the next flag
            } else if ((strcmp(argv[i], "-l" ) == 0) && (input_type == "")) {
                input_type = "-l";
                input_file = argv[i+1];
                single_file = false;
                i++;
            } else if ((strcmp(argv[i], "-o") == 0)) {
                output_file = argv[i+1];
                i++;
            } else if (strcmp(argv[i], "--peak_resolution") == 0) {
                TOLERANCE = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "-t") == 0) {
                THRESHOLD = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "--min_peakmatches") == 0) {
                MINMATCHES = stod(argv[i+1]);
                i++; 
            } else if (strcmp(argv[i], "--shifted_offset") == 0){
                shifted_offset = stoi(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "--edge_limit") == 0) {
                TOPProducts = stoi(argv[i+1]);
                i++; 
            } else if (strcmp(argv[i], "--cluster_info_in") == 0 ){
                cluster_info_in = argv[i+1];
                i++;
            } else if (strcmp(argv[i], "--connected_components_out") == 0 ){
                connected_components_out = argv[i+1];
                i++;
            } else if (i + 2 < argc) {
                if (strcmp(argv[i], "-f") == 0) {
                    SELECTIONRANGE = stod(argv[i+1]);
                    FILTERK = stod(argv[i+2]); 
                    i++;
                    i++;
                } 
            }
        } 
    }
    auto start = chrono::steady_clock::now();
    if (single_file) {
        vector<string> filenames;
        filenames.push_back(input_file);
        bool file_exist = examine_filenames(filenames);
        if (!file_exist) {
            exit(1);
        }
        pair<pair<table_t, table_t>, vector<PepMetadata>> input = parse_inputs_from_files(filenames);
        cout << "size: " << input.second.size() << endl;
        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        cout << "Input parsed: " << elapsed_seconds.count() << endl;
        start = chrono::steady_clock::now();
        unordered_map<int, unordered_map<int, double>> matches = all_versus_all(input.first, input.second, output_file, cluster_info_in, connected_components_out);
        end = chrono::steady_clock::now();
        elapsed_seconds = end-start;
        cout << "Processed data: " << elapsed_seconds.count() << endl;
    } else {
        vector<string> filenames = get_files(input_file);
        bool file_exist = examine_filenames(filenames);
        if (!file_exist) {
            exit(1);
        }
        pair<pair<table_t, table_t>, vector<PepMetadata>> input = parse_inputs_from_files(filenames);
        cout << "size: " << input.second.size() << endl;
        auto end = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        cout << "Input parsed: " << elapsed_seconds.count() << endl;
        start = chrono::steady_clock::now();
        unordered_map<int, unordered_map<int, double>> matches = all_versus_all(input.first, input.second, output_file, cluster_info_in, connected_components_out);
        end = chrono::steady_clock::now();
        elapsed_seconds = end-start;
        cout << "Processed data: " << elapsed_seconds.count() << endl;
    }
    return 0;
}
