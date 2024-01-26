#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include "bits/stdc++.h"
using namespace std;

#define FOR(i, hi) for (int i = 0; i < hi; ++i)
#define FORR(i, hi) for (int i = hi - 1; i >= 0; --i)
#define For(i, lo, hi) for (int i = lo; i < hi; ++i)
#define Forr(i, hi, lo) for (int i = hi - 1; i >= lo; --i)


typedef struct {
    //index of peak in the spectrum it came from
    int peak_idx;
    int scan;
    int local_cluster_idx;
    double pepmass;
    double true_mz;
    double mass;
} SpectrumInfo;

typedef struct {
    int mz_bin;
    int mass_bin;
    int scan;
    double true_mz;
    double mass;
} PeakInfo;

typedef struct {
    double pepmass;
    double RT;
    long double magnitude;
} Spec_info ;

typedef struct {
    int a, b;
    double massprod;
} Match;

typedef struct {
    int scan;
    double pepmass;
    double mz;
} top_four_Info;

typedef struct {
    int ind;
    double mass;
    double mz;
} Top_peak_t;

/*peakstorage datastructure**************************************************************/

typedef vector<vector<SpectrumInfo>> mztospectra_t;
/*index: the mz value bin index; value: the peaks with a mz value within the threshold*/

typedef unordered_map<int, vector<PeakInfo>> spectratomz_t;
/*index: the scan-number of each spectra; value: the peaks of that scan in the bin*/

typedef unordered_map<int, Spec_info> spectra_general_info_t;
/*index: the scan-number of each spectra; value: the general peak information (pepmass, RT) of that scan in the bin*/

typedef pair<spectratomz_t, spectra_general_info_t> spectratoinfo_t; 
/*the combination of peak informatin and general information for each scan*/


/*topfour datastructure ******************************************************************/

typedef unordered_map<int, vector<Top_peak_t>> scan_to_top_peaks_t;
/*index int: the scan-number of each spectra; values: a vector containing the top-four peaks of the scan*/

//typedef vector<top_four_Info> mz_to_topfour_t;
typedef vector<vector<int>> mz_to_topfour_t;
/*index int: the mz value of the topfour peaks; values: a vector containing the scans with at least one top-four peak in that bin*/

typedef vector<mz_to_topfour_t> pepmass_spectra_t;
/*index int: pepmass bins; values: the mz_to_topfour peaks vector with index as the mz bin value and the information stored as top-four peaks*/

typedef pair<pepmass_spectra_t, scan_to_top_peaks_t> top_four_combined_t;


/*cluster storage results*****************************************************************/
typedef unordered_map<int, int> scan_to_clustercenter_t;
/*index int: the cluster-center int; value: a vector containing the scans in that cluster; the first scan number refers to the cluster-center scan*/

typedef vector<vector<int>> clustercenter_to_scans_t;
/*index int: the scan number of spectras; value: the cluster-center of the scan*/

typedef pair<scan_to_clustercenter_t, clustercenter_to_scans_t> cluster_info_t;
/*the combination of the above datastructures*/

typedef pair<spectratoinfo_t, top_four_combined_t> topfour_pepmass_raw_t;

typedef vector<vector<int>> pepmass_distribution_t;

typedef pair<unordered_map<string, int>, unordered_map<int, string>> file_info_t;

static double TOLERANCE = 0.01;
static double THRESHOLD = 0.7;
static double MASSTOLERANCE = 1.0;
static double TOPFOURTOLERANCE = 0.1;
static int THRESHOLD_TOPFOUR = 5000;
static int THRESHOLD_INDEXING = 10;
static bool INDEXING_CHOICE = true;
static bool TOPFOUR_CHOICE = true;
static double SELECTIONRANGE = 50.0;
static int Output_cluster_minsize = 2;
static int FILTERK = 5;

vector<string> get_files(string content_files) {
    //filename should be name of file containing a path to an mgf file per line
    vector<string> filenames;
    ifstream f(content_files);
    string path_string;
    int load_tracer = 0;
    while ((f >> path_string)) {
        int path_len = path_string.length();
        int ext_begin = path_string.find_last_of(".");
        string file_ext = path_string.substr(ext_begin + 1);
        if ((file_ext == "mgf") || (file_ext == "mzML")) {
            filenames.push_back(path_string);
            load_tracer += 1;
        }
        /*
        if (path_len > 4) {
            string file_type = path_string.substr(path_len - 4, path_len);
            if (file_type.compare(".mgf") == 0) {
                filenames.push_back(path_string);
                load_tracer += 1;
            }
        }
        */
    }
    return filenames;
}


bool sort_tuple(const tuple<int, double, double>& a, const tuple<int, double, double>& b) {
    return (get<2>(a) > get<2>(b));
}

pair<file_info_t, pair<pepmass_distribution_t, topfour_pepmass_raw_t>> parse_inputs(vector<string> &filenames) {
    pair<file_info_t, pair<pepmass_distribution_t, topfour_pepmass_raw_t>> res;

    file_info_t &file_info = res.first;
    pepmass_distribution_t &pepmass_distribution = res.second.first;
    topfour_pepmass_raw_t &topfour_pepmass_raw = res.second.second;
    //pair<spectratoinfo_t, top_four_combined_t> res;

    unordered_map<string, int> &file_start_scan = file_info.first;
    unordered_map<int, string> &scan_file_src = file_info.second;

    spectratoinfo_t &spectratoinfo_shared = topfour_pepmass_raw.first;
    spectratomz_t &spectratomz_shared = spectratoinfo_shared.first;
    spectra_general_info_t &spectra_general = spectratoinfo_shared.second;

    top_four_combined_t &top_four_combined = topfour_pepmass_raw.second;
    pepmass_spectra_t &pepmass_spectra = top_four_combined.first;
    scan_to_top_peaks_t &scan_to_top_peaks = top_four_combined.second;

    pepmass_distribution.resize(500);
    pepmass_spectra.resize(500);

    //cout << "starting parsing inputs" << endl;
    int scan_tracer = 0;
    //start counting the files needed
    int verification = 0;
    for (int i=0; i<filenames.size(); i++) {
        string filename = filenames[i];
        int ext_begin = filename.find_last_of(".");
        string file_ext = filename.substr(ext_begin + 1);

        file_start_scan[filename] = scan_tracer;
        if (file_ext == "mgf") {
            ifstream f(filename);
            string a;
            while (getline(f, a)) {
                if (a == "BEGIN IONS") {
                    verification += 1;
                    double pepmass = -1;
                    double rtinsec = -1; 
                    while (getline(f, a)) {
                        if (a.rfind("PEPMASS", 0) == 0) {
                            pepmass = stod(a.substr(8));
                        }
                        if (a.rfind("RTINSECONDS", 0) == 0) {
                            rtinsec = stod(a.substr(12));
                        }
                        if ((pepmass != -1) && (rtinsec != -1)) {
                            /*
                            if (pepmass == -1) {
                                break;
                            }
                            */
                            vector<vector<tuple<int, double, double>>> curv_shared;
                            curv_shared.resize(20);
                            long double magnitude = 0;
                            bool nonempty = false;
                            /*
                            while (getline(f, a)) {
                                if (a == "END IONS") {
                                    break;
                                }
                                double mz, mass;
                                stringstream ss(a);
                                ss >> mz >> mass;
                                //make sure no two peaks are that close to each other 
                                if (abs(mz - prev_peak_mz) >= 2 * TOLERANCE) {
                                    // get bucket index for mz in both shared peaks settings
                                    int ind_shared = round(mz / TOLERANCE);
                                    tuple<int, double, double> elem{ind_shared, mz, sqrt(mass)};
                                    curv_shared.push_back(elem);
                                    magnitude += mass;
                                    prev_peak_mz = mz;
                                }
                                nonempty = true;
                            }
                            */
                            while (getline(f, a)) {
                                if (a == "END IONS") {
                                    break;
                                }
                                double mz, mass;
                                stringstream ss(a);
                                ss >> mz >> mass;
                                // get bucket index for mz in both shared peaks settings
                                int ind_shared = round(mz / TOLERANCE);

                                //Modification
                                int filter_k_index = round(mz / SELECTIONRANGE);
                                if (filter_k_index >= curv_shared.size()) {
                                    curv_shared.resize(filter_k_index * 3 / 2 + 1);
                                }
                                tuple<int, double, double> elem{ind_shared, mz, mass};
                                curv_shared[filter_k_index].push_back(elem);
                                //magnitude += mass;
                                nonempty = true;
                            }
                            vector<tuple<int, double, double>> filtered_shared_peaks;
                            //MODIFICATION FOR FILTERING TOPK
                            for (int filteringk_index=0; filteringk_index < curv_shared.size(); filteringk_index ++) {
                                vector<tuple<int, double, double>> primitive_bin = curv_shared[filteringk_index];
                                if (primitive_bin.size() >= FILTERK) {
                                    sort(primitive_bin.begin(), primitive_bin.end(), sort_tuple);
                                    curv_shared[filteringk_index] = {primitive_bin.begin(), primitive_bin.begin()+FILTERK};
                                }
                                int peaks_in_bin = primitive_bin.size();
                                int peak_preserved = min(peaks_in_bin, FILTERK);
                                for (int j = 0; j < peak_preserved; j++) {
                                    filtered_shared_peaks.push_back(primitive_bin[j]);
                                    double mass = get<2>(primitive_bin[j]);
                                    magnitude += mass * mass;
                                }
                            }
                            //cout << "start adding peaks" << endl;
                            if (nonempty) {
                            	vector<Top_peak_t> top_four_peaks;
                            	for (int j=0; j<4; j++) {
                            		top_four_peaks.push_back({0, 0.0, 0.0});
                            	}
                                //cout << "top_four initialized" << endl;
                                int scan = scan_tracer;
                                scan_tracer += 1;
                                magnitude = sqrt(magnitude);
                                int pepmass_idx = round(pepmass / MASSTOLERANCE); 
                                for (int k = 0; k < filtered_shared_peaks.size(); ++k) {
                                    tuple<int, double, double> p = filtered_shared_peaks[k];
                                    int ind = get<0>(p);
                                    double mz = get<1>(p);
                                    double mass = get<2>(p);
                                    //cout << "updating spectratomz once" << endl;
                                    spectratomz_shared[scan].push_back({ind, pepmass_idx, scan, mz, (double)(mass / magnitude)});
                                    int top_four_ind = round(mz / TOPFOURTOLERANCE);
                                    if (mass > top_four_peaks[2].mass) {
                                    	if (mass > top_four_peaks[1].mass) {
                                    		if (mass > top_four_peaks[0].mass) {
                                    			top_four_peaks[3] = {top_four_peaks[2].ind, top_four_peaks[2].mass, top_four_peaks[2].mz};
                                    			top_four_peaks[2] = {top_four_peaks[1].ind, top_four_peaks[1].mass, top_four_peaks[1].mz};
                                    			top_four_peaks[1] = {top_four_peaks[0].ind, top_four_peaks[0].mass, top_four_peaks[0].mz};
                                    			top_four_peaks[0] = {top_four_ind, mass, mz};
                                    		} else {
                                    			top_four_peaks[3] = {top_four_peaks[2].ind, top_four_peaks[2].mass, top_four_peaks[2].mz};
                                    			top_four_peaks[2] = {top_four_peaks[1].ind, top_four_peaks[1].mass, top_four_peaks[1].mz};
                                    			top_four_peaks[1] = {top_four_ind, mass, mz};
                                    		}

                                    	} else {
                                    		top_four_peaks[3] = {top_four_peaks[2].ind, top_four_peaks[2].mass, top_four_peaks[2].mz};
                                    		top_four_peaks[2] = {top_four_ind, mass, mz};
                                    	}

                                    } else {
                                    	if (mass > top_four_peaks[3].mass) {
                                    		top_four_peaks[3] = {top_four_ind, mass, mz};
                                    	}
                                    }
                                }
                                if (pepmass_idx >= pepmass_spectra.size()){
                                    pepmass_spectra.resize(pepmass_idx * 3 / 2 + 1);
                                    pepmass_distribution.resize(pepmass_idx * 3 / 2 + 1);
                                }

                                int current_top_four_bins = pepmass_spectra[pepmass_idx].size();
                                //add top-peaks to data-structure
                                //cout << "finish adding peaks" << endl;
                                for (int j=0; j<4; j++) {
                               		int top_four_ind = top_four_peaks[j].ind;
                                    double top_four_intensity = top_four_peaks[j].mass;
                                    double top_four_mz = top_four_peaks[j].mz;
                                    if (top_four_intensity > 0.0){
                                        if (top_four_ind >= current_top_four_bins){
                                            pepmass_spectra[pepmass_idx].resize(top_four_ind * 3 / 2 + 1);
                                        }
                                        //pepmass_spectra[pepmass_idx][top_four_ind].push_back({scan, pepmass, top_four_mz});
                                        pepmass_spectra[pepmass_idx][top_four_ind].push_back(scan);
                                        scan_to_top_peaks[scan] = top_four_peaks;
                                    }

                            	}
                                scan_to_top_peaks[scan] = top_four_peaks;
                                pepmass_distribution[pepmass_idx].push_back(scan);
                                spectra_general[scan] = {pepmass, rtinsec, magnitude};
                                scan_file_src[scan] = filename;
                                //cout << "a spectra parsed" << endl;
                            }
                            break;
                        }  
                    }
                }
            }
            f.close();
        }
    }
    return res;
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

bool cmp(Match x, Match y) {
    return (x.massprod > y.massprod);
}

void update_product(mztospectra_t &local_mztospectra_center, vector<double> &product_res,  vector<PeakInfo> &peaks_share, double pepmass_scan) {
    int pepmass_idx = round(pepmass_scan / MASSTOLERANCE);
    for (int j = 0; j < peaks_share.size(); ++j) {
        PeakInfo peak_shared = peaks_share[j];
        int bin_shared = peak_shared.mz_bin;
        double mz_shared = peak_shared.true_mz;
        double mass = peak_shared.mass;
        bool used = false;
        if (bin_shared < local_mztospectra_center.size()) {
            for(auto t1 = local_mztospectra_center[bin_shared].begin(); t1 != local_mztospectra_center[bin_shared].end(); ++t1) {
                if (abs(pepmass_scan - t1 -> pepmass) <= MASSTOLERANCE) {
                    product_res[t1->local_cluster_idx] += mass * t1->mass;
                }
            }
        }
        if ((bin_shared > 0) && (bin_shared - 1 < local_mztospectra_center.size())) {
            for(auto t2 = local_mztospectra_center[bin_shared - 1].begin(); t2 != local_mztospectra_center[bin_shared - 1].end(); ++t2) {
                if (abs(pepmass_scan - t2 -> pepmass) <= MASSTOLERANCE) {
                    if(abs(mz_shared - t2->true_mz) <= TOLERANCE) {
                        product_res[t2->local_cluster_idx] += mass * t2->mass;
                    }
                }
            }
        }
        if (bin_shared + 1 < local_mztospectra_center.size()) {
            for(auto t3 = local_mztospectra_center[bin_shared + 1].begin(); t3 != local_mztospectra_center[bin_shared + 1].end(); ++t3) {
                if (abs(pepmass_scan - t3 -> pepmass) <= MASSTOLERANCE) {
                    if(abs(mz_shared - t3->true_mz) <= TOLERANCE) {
                        product_res[t3->local_cluster_idx] += mass * t3->mass;
                    }
                }
            }
        }
    }
}



int brute_force_clustering(spectratoinfo_t& spectrainfo_all, vector<int>& scan_numbers_in_bin, cluster_info_t& cluster_info) {
    if (scan_numbers_in_bin.size() == 0) {
        return 0;
    }
    int total_pairwise_counter = 0;
    spectratomz_t &spectratomz_shared = spectrainfo_all.first;
    spectra_general_info_t &spectra_general = spectrainfo_all.second;

    scan_to_clustercenter_t &spectra_cluster = cluster_info.first;
    clustercenter_to_scans_t &cluster_content = cluster_info.second;

    unordered_map<int, bool> cluster_marker;
    spectratomz_t local_spectratomz_center;
    for (int scan_idx = 0; scan_idx < scan_numbers_in_bin.size(); scan_idx ++) {
        int fixed_scan = scan_numbers_in_bin[scan_idx];
        vector<PeakInfo>& peaks_share = spectratomz_shared[fixed_scan];
        if (spectra_cluster[fixed_scan] >= 0) {
            int cluster_idx_global = spectra_cluster[fixed_scan];
            int cluster_center_scan = cluster_content[cluster_idx_global][0];
            if (!cluster_marker[cluster_idx_global]) {
                local_spectratomz_center[cluster_idx_global] = spectratomz_shared[cluster_center_scan];
                cluster_marker[cluster_idx_global] = true;
            }
        } else {
            int cluster_belonging = -1;
            double max_prod = 0.0;
            for (auto iter = local_spectratomz_center.begin(); iter != local_spectratomz_center.end(); iter ++) {
                vector<Match> matches;
                int cluster_idx_global = iter -> first;
                vector<PeakInfo>& cluster_peaks = iter -> second;
                int cluster_center_scan = cluster_content[cluster_idx_global][0];
                if (abs(spectra_general[cluster_center_scan].pepmass - spectra_general[fixed_scan].pepmass) > MASSTOLERANCE) {
                    continue;
                }
                double product_res = 0.0;
                for (int i = 0; i < peaks_share.size(); i++) {
                    PeakInfo peak_shared = peaks_share[i];
                    for (int j = 0; j < cluster_peaks.size(); j++) {
                        PeakInfo peak_cluster = cluster_peaks[j];
                        if (abs(peak_shared.true_mz - peak_cluster.true_mz) <= TOLERANCE) {
                            product_res += peak_shared.mass * peak_shared.mass;
                        } 
                    }
                }
                if (product_res > max_prod) {
                    max_prod = product_res;
                    cluster_belonging = cluster_idx_global;
                }
            }

            if ((max_prod > 0.0) && (cluster_belonging >= 0)) {
                cluster_content[cluster_belonging].push_back(fixed_scan);
                spectra_cluster[fixed_scan] = cluster_belonging;
                cluster_marker[cluster_belonging] = true;
            } else {
                //the cluster is brand new , creating a new cluster
                int cluster_idx_global = cluster_content.size();
                //int local_cluster_idx = query_specs.size();
                //query_specs.push_back(cluster_idx_global);
                vector<int> new_cluster_sublist;
                cluster_content.push_back(new_cluster_sublist);
                cluster_content[cluster_idx_global].push_back(fixed_scan);
                spectra_cluster[fixed_scan] = cluster_idx_global;
                cluster_marker[cluster_idx_global] = true;
                local_spectratomz_center[cluster_idx_global] = peaks_share;
            }
        }
    }
    return total_pairwise_counter;
}



int cluster_local_bin(spectratoinfo_t& spectrainfo_all, vector<int>& scan_numbers_in_bin, cluster_info_t& cluster_info) {
    //the spectratomz of all spectras
    if (scan_numbers_in_bin.size() == 0) {
        return 0;
    }
    int total_pairwise_counter = 0;
    spectratomz_t &spectratomz_shared= spectrainfo_all.first;

    //the general information of all spectras
    spectra_general_info_t &spectra_general= spectrainfo_all.second;

    scan_to_clustercenter_t& spectra_cluster = cluster_info.first;

    clustercenter_to_scans_t& cluster_content = cluster_info.second;

    mztospectra_t local_mztospectra_center;
    unordered_map<int, bool> cluster_marker;
    spectratomz_t local_spectratomz_center;
    vector<int> query_specs;

    local_mztospectra_center.resize(200000);

    for (int scan_idx = 0; scan_idx < scan_numbers_in_bin.size(); scan_idx ++) {
        int fixed_scan = scan_numbers_in_bin[scan_idx];
        vector<PeakInfo>& peaks_share = spectratomz_shared[fixed_scan];
        double pepmass_scan = spectra_general[fixed_scan].pepmass;

        if (spectra_cluster[fixed_scan] >= 0) {
            int cluster_idx_global = spectra_cluster[fixed_scan]; 
            int cluster_center_scan = cluster_content[cluster_idx_global][0];
            if (!cluster_marker[cluster_idx_global]) {
                //add local_mztospectra_center 
                int peak_idx = 0;

                local_spectratomz_center[cluster_idx_global] = spectratomz_shared[cluster_center_scan];
                cluster_marker[cluster_idx_global] = true;
                int local_cluster_idx = query_specs.size();
                query_specs.push_back(cluster_idx_global);

                for (auto iter = peaks_share.begin(); iter != peaks_share.end(); ++iter) {
                    int ind_peak = iter->mz_bin;
                    int scan = iter->scan;
                    double mass = iter->mass;
                    double mz = iter->true_mz;
                    if (ind_peak >= local_mztospectra_center.size()) {
                        local_mztospectra_center.resize(ind_peak * 3 / 2 + 1);
                    }
                    local_mztospectra_center[ind_peak].push_back({peak_idx, cluster_center_scan, local_cluster_idx, pepmass_scan, mz, mass});
                    peak_idx += 1;
                }
            }
        } else {
            vector<double> product_res;
            product_res.resize(local_spectratomz_center.size());
            double pepmass_scan = spectra_general[fixed_scan].pepmass;
            update_product(local_mztospectra_center, product_res, peaks_share, pepmass_scan);
            int cluster_belonging = -1;
            double max_prod = 0.0;

            for (int local_cluster_idx = 0; local_cluster_idx < product_res.size(); local_cluster_idx++) {
                double prod_tmp = product_res[local_cluster_idx];
                if ((prod_tmp >= THRESHOLD) && (prod_tmp > max_prod)) {
                    max_prod = prod_tmp;
                    int cluster_idx_global = query_specs[local_cluster_idx];
                    cluster_belonging = cluster_idx_global;
                } 
            }
            if ((max_prod > 0.0) && (cluster_belonging >= 0)) {
                cluster_content[cluster_belonging].push_back(fixed_scan);
                spectra_cluster[fixed_scan] = cluster_belonging;
                cluster_marker[cluster_belonging] = true;
            } else {
                //the cluster is brand new , creating a new cluster
                int cluster_idx_global = cluster_content.size();
                int local_cluster_idx = query_specs.size();
                query_specs.push_back(cluster_idx_global);
                vector<int> new_cluster_sublist;
                cluster_content.push_back(new_cluster_sublist);
                cluster_content[cluster_idx_global].push_back(fixed_scan);
                spectra_cluster[fixed_scan] = cluster_idx_global;
                cluster_marker[cluster_idx_global] = true;


                local_spectratomz_center[cluster_idx_global] = peaks_share;
                
                int peak_idx = 0;
                for (auto iter = peaks_share.begin(); iter != peaks_share.end(); ++iter) {
                    int ind_peak = iter->mz_bin;
                    int scan = iter->scan;
                    double mass = iter->mass;
                    double mz = iter->true_mz;
                    if (ind_peak >= local_mztospectra_center.size()) {
                        local_mztospectra_center.resize(ind_peak * 3 / 2 + 1);
                    }
                    local_mztospectra_center[ind_peak].push_back({peak_idx, fixed_scan, local_cluster_idx, pepmass_scan, mz, mass});
                    peak_idx += 1;
                }
            }
        }
    }
    return total_pairwise_counter;
}




void generate_clusters(pepmass_distribution_t &pepmass_distribution, spectratoinfo_t& spectrainfo_all, top_four_combined_t& top_four_combined, cluster_info_t& cluster_info) {
    pepmass_spectra_t &pepmass_spectra = top_four_combined.first;
    scan_to_top_peaks_t &scan_to_top_peaks = top_four_combined.second;

    scan_to_clustercenter_t& spectra_cluster = cluster_info.first;
    clustercenter_to_scans_t& cluster_content = cluster_info.second;

    spectratomz_t &spectratomz_all = spectrainfo_all.first;
    for (auto it = spectratomz_all.begin(); it != spectratomz_all.end(); ++it) {
        int query_scan = it->first;
        spectra_cluster[query_scan] = -1;
    }

    float total_time = 0.0;
    float pepmassbin_abovethresh_total = 0.0; 
    int bins_with_topfour = 0;
    //long total_pairwise_prod_topfour = 0;

    cout << "total number of pepmass bins: " << pepmass_spectra.size() << endl;
    for (int pepmass_bin = 0; pepmass_bin < pepmass_spectra.size(); pepmass_bin ++) {
        //cout << "start on pepmass_bin " << pepmass_bin << endl;
        int spectra_counter_for_pepmass = pepmass_distribution[pepmass_bin].size(); 
        if (spectra_counter_for_pepmass > THRESHOLD_TOPFOUR) {
            //case 1: the total number of spectra in the pepmass bin is worthy to use topfour
            //int pairwise_prod_pepmass_bin = 0;
            //cout << "entering case 1 " << endl;
            auto topfour_start = chrono::steady_clock::now();
            mz_to_topfour_t& mz_to_topfour = pepmass_spectra[pepmass_bin];
            for (int top_four_bin = 0; top_four_bin < mz_to_topfour.size(); top_four_bin ++) {
                vector<int> scan_numbers_in_topfour_bin; 
                if (pepmass_bin > 0) {
                    mz_to_topfour_t& mz_to_topfour_prev = pepmass_spectra[pepmass_bin - 1];
                    if (top_four_bin < mz_to_topfour_prev.size()){
                        scan_numbers_in_topfour_bin.insert(scan_numbers_in_topfour_bin.end(), mz_to_topfour_prev[top_four_bin].begin(), mz_to_topfour_prev[top_four_bin].end());
                    }
                } 
                scan_numbers_in_topfour_bin.insert(scan_numbers_in_topfour_bin.end(), mz_to_topfour[top_four_bin].begin(), mz_to_topfour[top_four_bin].end());
                //brute_force_clustering(spectrainfo_all, scan_numbers_in_topfour_bin, cluster_info);
                if ((INDEXING_CHOICE == true) && (scan_numbers_in_topfour_bin.size() > THRESHOLD_INDEXING)) {
                    cluster_local_bin(spectrainfo_all, scan_numbers_in_topfour_bin, cluster_info);
                } else {
                    brute_force_clustering(spectrainfo_all, scan_numbers_in_topfour_bin, cluster_info);
                }
                
            }
            auto topfour_end = chrono::steady_clock::now();
            chrono::duration<double> topfour_pepmass_bin_time = topfour_end - topfour_start;
            pepmassbin_abovethresh_total += float(topfour_pepmass_bin_time.count());
            bins_with_topfour += 1;
        } else {
            vector<int> scan_numbers_in_pepmassbin;
            if (pepmass_bin > 0) {
                vector<int> &prev_scan_numbers_in_pepmassbin = pepmass_distribution[pepmass_bin - 1];
                scan_numbers_in_pepmassbin.insert(scan_numbers_in_pepmassbin.end(), prev_scan_numbers_in_pepmassbin.begin(), prev_scan_numbers_in_pepmassbin.end());
            }
            vector<int> &current_scan_numbers_in_pepmassbin = pepmass_distribution[pepmass_bin];
            scan_numbers_in_pepmassbin.insert(scan_numbers_in_pepmassbin.end(), current_scan_numbers_in_pepmassbin.begin(), current_scan_numbers_in_pepmassbin.end());

            if ((INDEXING_CHOICE == true) && (scan_numbers_in_pepmassbin.size() > THRESHOLD_INDEXING)) {
                brute_force_clustering(spectrainfo_all, scan_numbers_in_pepmassbin, cluster_info);
            } else {
                cluster_local_bin(spectrainfo_all, scan_numbers_in_pepmassbin, cluster_info);
            }
        }
    }
    cout << "the total time needed for top four methods on pepmass bins above the threshold: " << pepmassbin_abovethresh_total << endl;
    cout << "the total number of pepmass bins with topfour applied: "<<  bins_with_topfour << endl;
}

void write_cluster_centers(spectratoinfo_t& spectrainfo_all, vector<int>& scans, string cluster_center_res) {
    spectratomz_t &spectratomz_shared= spectrainfo_all.first;
    spectra_general_info_t &spectra_general_info= spectrainfo_all.second;
    ofstream cluster_center_writer;
    cluster_center_writer.open(cluster_center_res);
    for (int i=0; i < scans.size(); i++) {
        int center_scan = scans[i];
        vector<PeakInfo>& cluster_peaks = spectratomz_shared[center_scan];
        double pepmass = spectra_general_info[center_scan].pepmass;
        double RT = spectra_general_info[center_scan].RT;
        long double magnitude = spectra_general_info[center_scan].magnitude;
        cluster_center_writer << "BEGIN IONS"<< endl;
        cluster_center_writer << "PEPMASS=" << pepmass << endl;
        cluster_center_writer << "RTINSECONDS=" << RT << endl;
        for (int j = 0; j < cluster_peaks.size(); j++) {
            double true_mz = cluster_peaks[j].true_mz;
            double mass = cluster_peaks[j].mass;
            double true_mass = (mass * magnitude) * (mass * magnitude);
            cluster_center_writer << true_mz << " " << true_mass << endl;
        }
        cluster_center_writer << "END IONS" << endl;

    }
    cluster_center_writer.close();
}


int main(int argc, char* argv[]) {
    string input_type = "";
    string output_type = "";
    string input_file = "";
    string output_file = "";
    string arguments = "";
    string cluster_info_output = "cluster_info.csv";
    string printed_cluster_centers = "centers.mgf"; 
    for (int k = 0; k < argc; k++) {
        arguments.append(argv[k]);
        arguments.append(" ");
    }

    cout << arguments << endl;
    bool single_file = true;
    if (argc < 5) {
        cout << "instructions for running the code: " << endl;
        cout << "usage: ./clustering [-i/-l <input_path>] -o <output_path> [--no_topfour] [--bruteforce] [--topfour_threshold <default=5000>] --topfour_resolution <default=0.1> --pepmass_resolution <default=1.0> --peak_resolution <default=0.01> --product_threshold <default = 0.7>" << endl;
        cout << "\n"<< endl;
        cout << "--help: print instructions"<< endl;
        cout << "-i <input_path> : cluster spectra from a single input file " << endl; 
        cout << "-l <input_path> : cluster spectra from a list of input files" << endl; 
        cout << "-o <output_path>: output path for clustering results" << endl;
        cout << "-t <default=0.7> (optional): the minimum threshold of dot product between a cluster member and its corresponding cluster center"<< endl;
        cout << "-c <default=centers.mgf> (optional): the cluster centers produced by the code" << endl;
        cout << "-f X<default=50.0> Y<default=5> (optional): add filtering to peaks in the input spectra. X: the range of filtering peaks; Y: preserving the top Y peaks in each mz range X during filtering" << endl;
        cout << "-s <default=cluster_info.csv>: the information each cluster contains" << endl;
        cout << "--no_topfour (optional): don't apply top four peaks in clustering" << endl;
        cout << "--bruteforce (optional): use bruteforce technique rather than indexing in clustering" << endl;
        cout << "--topfour_threshold <default=5000> (optional): the minimum number of spectra in a pepmass bin to trigger top four peaks filtering. Can't coexist with flag --no_topfour" << endl;
        cout << "--topfour_resolution <default=0.1> (optional): the mz tolerance of two top four peaks. Can't coexist with flag --no_topfour" << endl;
        cout << "--pepmass_resolution <default=1.0> (optional): the minimum pepmass difference to perform dot product between two spectra"<< endl;
        cout << "--peak_resolution <default=0.01> (optional): the minimum mz difference for adding the product between two peaks in calculating the spectra pairwise dotproduct"<< endl;
        cout << "--cluster_minsize <default=2>: the minimum cluster size for a cluster center to be preserved for spectra networking"<< endl; 
        cout << "\n" << endl;
        return 0; 
    }
    for (int i = 1; i < argc; i++) {
        //cout << argc << " " << argv[i] << endl;
        if ((i == argc - 1) && (i==1)) {
            //cout << argc << " " << argv[i] << endl;
            if (strcmp(argv[i], "--help") == 0) {
                cout << "instructions for running the code: " << endl;
                cout << "usage: ./clustering [-i/-l <input_path>] -o <output_path> [--no_topfour] [--bruteforce] [--topfour_threshold <default=5000>] --topfour_resolution <default=0.1> --pepmass_resolution <default=1.0> --peak_resolution <default=0.01> --product_threshold <default = 0.7>" << endl;
                cout << "\n"<< endl;
                cout << "--help: print instructions"<< endl;
                cout << "-i <input_path> : cluster spectra from a single input file " << endl; 
                cout << "-l <input_path> : cluster spectra from a list of input files" << endl; 
                cout << "-o <output_path>: output path for clustering results" << endl;
                cout << "-t <default=0.7> (optional): the minimum threshold of dot product between a cluster member and its corresponding cluster center"<< endl;
                cout << "-c <default=centers.mgf> (optional): the cluster centers produced by the code" << endl;
                cout << "-f X<default=50.0> Y<default=5> (optional): add filtering to peaks in the input spectra. X: the range of filtering peaks; Y: preserving the top Y peaks in each mz range X during filtering" << endl;
                cout << "-s <default=cluster_info.csv>: the information each cluster contains" << endl; 
                cout << "--no_topfour (optional): don't apply top four peaks in clustering" << endl;
                cout << "--bruteforce (optional): use bruteforce technique rather than indexing in clustering" << endl;
                cout << "--topfour_threshold <default=5000> (optional): the minimum number of spectra in a pepmass bin to trigger top four peaks filtering. Can't coexist with flag --no_topfour" << endl;
                cout << "--topfour_resolution <default=0.1> (optional): the mz tolerance of two top four peaks. Can't coexist with flag --no_topfour" << endl;
                cout << "--pepmass_resolution <default=1.0> (optional): the minimum pepmass difference to perform dot product between two spectra"<< endl;
                cout << "--peak_resolution <default=0.01> (optional): the minimum mz difference for adding the product between two peaks in calculating the spectra pairwise dotproduct"<< endl;
                cout << "--cluster_minsize <default=2>: the minimum cluster size for a cluster center to be preserved for spectra networking"<< endl; 
                cout << "\n" << endl;
                return 0; 
            }
        }

        if (strcmp(argv[i], "--no_topfour") == 0) {
            TOPFOUR_CHOICE = false;
        } else if (strcmp(argv[i], "--bruteforce")) {
            INDEXING_CHOICE = false;
        }

        if (i + 1 != argc) {
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
            } else if ((strcmp(argv[i], "--topfour_threshold") == 0)) {
                if (TOPFOUR_CHOICE == false) {
                    cout << argv[i] <<" " << argv[i+1] <<": unable to use top-four-peak filtering" << endl;
                    cout << "flag --no_topfour can't coexist with flag --topfour_threshold" << endl;
                    return 1;
                }
                THRESHOLD_TOPFOUR = stoi(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "--topfour_resolution") == 0) {
                if (TOPFOUR_CHOICE == false) {
                    cout << argv[i] <<" " << argv[i+1] <<": unable to use top-four-peak filtering" << endl;
                    cout << "flag --no_topfour can't coexist with flag --topfour_threshold" << endl;
                    return 1;
                }
                TOPFOURTOLERANCE = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "--pepmass_resolution") == 0) {
                MASSTOLERANCE = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "--peak_resolution") == 0) {
                TOLERANCE = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "-t") == 0) {
                THRESHOLD = stod(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "-c") == 0) {
                printed_cluster_centers = argv[i+1];
                i++;
            } else if (strcmp(argv[i], "--cluster_minsize") == 0) {
                Output_cluster_minsize = stoi(argv[i+1]);
                i++;
            } else if (strcmp(argv[i], "-s") == 0) {
                cluster_info_output =  argv[i+1];
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

    if (input_file == "") {
        cout << "missing argument, use -i <input_file> or -l <input_file> to specify the path to spectras to cluster"<< endl;
        return 1;
    }

    if (output_file == "") {
        cout << "missing argument, use -o <output_file> to specify the path of the output_file"<< endl;
        return 1;
    }

    if (TOPFOUR_CHOICE == false) {
        THRESHOLD_TOPFOUR = INT_MAX;
    }

    /*
    if (argc < 5) {
        cout << "insufficient number of variables" << endl;
        cout << "correct format1: ./<executable> -i <input file> -o <output file>" << endl;
        cout << "correct format2: ./<executable> -l <input list> -o <output file>" << endl;
        exit(1);
    }

    if(argc > 1) {
        input_type = argv[1];
        if (input_type == "-i") {
            single_file = true;
        } else if (input_type == "-l") {
            single_file = false;
        } else {
            cout << input_type <<" unidentified, requires -i or -l flags in this location" << endl;
            cout << "correct format1: ./<executable> -i <input file> -o <output file>" << endl;
            cout << "correct format2: ./<executable> -l <input list> -o <output file>" << endl;
            exit(1);
        }
    }
    if(argc > 2) { 
        input_file = argv[2];
    }

    if(argc > 3) {
        output_type = argv[3];
        if (output_type != "-o") {
            cout << output_type <<" unidentified, requires -o flag in this location" << endl;
            cout << "correct format1: ./<executable> -i <input file> -o <output file>" << endl;
            cout << "correct format2: ./<executable> -l <input list> -o <output file>" << endl;
            exit(1);
        }
    }
    
    if(argc > 4) { 
        output_file = argv[4];
    }

    if (argc > 5) {
        string extra_info = argv[5];
        cout << extra_info <<" unidentified variable" << endl;
        cout << "correct format1: ./<executable> -i <input file> -o <output file>" << endl;
        cout << "correct format2: ./<executable> -l <input list> -o <output file>" << endl;
        exit(1);
    } 
    */

    vector<string> filenames;
    if (single_file) {
        filenames.push_back(input_file);
        cout << "input from single file " << endl;
    } else {
        filenames = get_files(input_file);
        cout << "input from file list " << endl;
    }
    
    
    vector<int> scan_centers;
    int parsed_file = 0;
    int start_scan = 0;
    int file_idx_tracer = 0;

    //initialize the cluster center
    auto cluster_start = chrono::steady_clock::now();
    //cout << "tag1" << endl;
    auto preprocessing_start = chrono::steady_clock::now();
    pair<file_info_t ,pair<pepmass_distribution_t, topfour_pepmass_raw_t>> parsed_input = parse_inputs(filenames);
    //cout << "tag2" << endl;

    auto preprocessing_end = chrono::steady_clock::now();
    chrono::duration<double> preprocessing_time = preprocessing_end - preprocessing_start;
    cout << "preprocessing_time: " << preprocessing_time.count() << endl;
    //cout << "tag3" << endl;

    auto generating_clusters_start = chrono::steady_clock::now();

    file_info_t &file_info = parsed_input.first;
    unordered_map<string, int> &file_start_scan = file_info.first;
    unordered_map<int, string> &scan_file_src = file_info.second;
    //cout << "tag4" << endl;    

    pepmass_distribution_t &pepmass_distribution = parsed_input.second.first;
    topfour_pepmass_raw_t &topfour_pepmass_raw = parsed_input.second.second;
    cluster_info_t cluster_info;
    spectratoinfo_t& spectrainfo_all = topfour_pepmass_raw.first;
    top_four_combined_t& top_four_combined = topfour_pepmass_raw.second;
    //cout << "tag5" << endl;
    generate_clusters(pepmass_distribution, spectrainfo_all, top_four_combined, cluster_info);
    //cout << "tag6" << endl;
    auto generating_clusters_end = chrono::steady_clock::now();
    chrono::duration<double> generating_clusters_time = generating_clusters_end - generating_clusters_start;
    cout << "time for generating clusters " << generating_clusters_time.count() << endl;
    cout << "total number of spectra: " <<  cluster_info.first.size() << endl;
    cout << "total number of clusters: " << cluster_info.second.size() << endl;


    auto cluster_end = chrono::steady_clock::now();
    chrono::duration<double> total_cluster_time = cluster_end - cluster_start;
    cout << "total clustering time: " << total_cluster_time.count() << endl;

    auto writing_begin = chrono::steady_clock::now();
    ofstream outf;
    outf.open(output_file);

    ofstream cluster_outputf;
    cluster_outputf.open(cluster_info_output);
    
    cluster_outputf << "cluster_idx" << "\t" << "average mz"<< "\t" << "average RT" << "\t"  << "num spectra"  << endl;
    outf << "cluster_idx" << "\t" << "scan"<< "\t" << "mz" << "\t" << "RTINSECONDS" << "\t" << "index in file" << "\t" <<"source filename" << endl;
    
    scan_to_clustercenter_t& spectra_cluster = cluster_info.first;
    clustercenter_to_scans_t& cluster_content = cluster_info.second;
    spectra_general_info_t &spectra_general = spectrainfo_all.second;
    for (int cluster_idx = 0; cluster_idx < cluster_content.size(); cluster_idx ++) {
        vector<int> content_tmp = cluster_content[cluster_idx];
        if (content_tmp.size() < Output_cluster_minsize) {
            continue;
        }
        double pepmass_sum = 0.0;
        double RT_sum = 0.0;
        double intensity_sum ;
        int num_spectra = content_tmp.size();
        int connected_component_index;

        for (int content_idx = 0; content_idx < content_tmp.size(); content_idx ++) {
            int query_scan = content_tmp[content_idx];
            if ((content_idx == 0) && (content_tmp.size() >= Output_cluster_minsize))  {
                scan_centers.push_back(query_scan);
            }
            double mz = spectra_general[query_scan].pepmass;
            double RT = spectra_general[query_scan].RT;

            pepmass_sum += mz;
            RT_sum += RT;
            
            string src_filename = scan_file_src[query_scan];
            int local_index = query_scan - file_start_scan[src_filename];

            outf << cluster_idx << "\t" << query_scan << "\t" << mz << "\t" << RT << "\t" << local_index << "\t" << src_filename << endl;
        }
        double pepmass_sum_avg = pepmass_sum / num_spectra;
        double RT_avg = RT_sum / num_spectra;
        cluster_outputf << cluster_idx << "\t" << pepmass_sum_avg << "\t" << RT_avg << "\t" << num_spectra <<  endl;

    }

    outf.close();
    cluster_outputf.close();

    write_cluster_centers(spectrainfo_all, scan_centers, printed_cluster_centers);
    auto writing_end = chrono::steady_clock::now();
    chrono::duration<double> writing_duration = writing_end - writing_begin;
    cout << "Output writing time: " << writing_duration.count() << endl;
    return 0;
}
