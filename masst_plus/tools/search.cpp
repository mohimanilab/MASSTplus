#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <utility>
#include <fstream>
#include <experimental/filesystem>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../mgf/spectra_mgf.hpp"
#include "../mgf/parse_mgf.hpp"
#include "../spectral_index.hpp"

namespace fs = std::experimental::filesystem;

void print_help_message(){
	printf("Usage: ./search <query_file(list)_name>\n");
	printf("       -h: print this help message\n");
	printf("       --query-list, -q: indicate query is a file list\n");
	printf("       -t <threshold>: specify threshold for search, default 0.8\n");
  printf("       -a <analog>: search all precursor masses with and without mass shift\n");
  printf("       -p <peak-tolerance>: peak tolerance, default 0.02\n");
	return;
}

void parse_args(int argc, char** argv,
				std::string& query_file_name,
				bool& query_is_list,
				double& threshold,
        bool& error_tolerant,
				bool& hybrid,
				double& precursor_tolerance,
        double& peak_tolerance,
				std::string& library_path,
				std::string& output_path,
				bool& debug){
	static struct option long_options[] = {
		{"reference-list", no_argument, 0, 'r'},
    {"thresh", required_argument, 0, 't'},
    {"analog", no_argument, 0, 'a'},
		{"hybrid", no_argument, 0, 'x'},
		{"debug", no_argument, 0, 'd'},
    {"peaktol", required_argument, 0, 'p'},
		{"precursortol", required_argument, 0, 'q'},
		{"library", required_argument, 0, 'l'},
		{"output", required_argument, 0, 'o'}
  };
	char c;
	int option_index = 0;

	while((c = getopt_long(argc, argv, "raxdt:p:q:l:o:h", long_options, &option_index)) != -1){
		switch(c){
			case 'r':
				query_is_list = true;
				break;
			case 't':
				threshold = atof(optarg);
        break;
      case 'a':
				error_tolerant = true;
        break;
			case 'x':
				hybrid = true;
        break;
      case 'p':
        peak_tolerance = atof(optarg);
        break;
			case 'q':
        precursor_tolerance = atof(optarg);
        break;
			case 'l':
				library_path = optarg;
				break;
			case 'o':
				output_path = optarg;
				break;
			case 'd':
				debug = true;
				break;
			case 'h':
				print_help_message();
				exit(0);
			default:
				abort();
		}
	}

	if (argc - optind != 1){
		printf("False number of arguments!\n");
		print_help_message();
		return;
	}
	query_file_name = std::string(argv[optind++]);

	return;
}

// find the "best" set of matches (this is a heuristic)
void getTopMatches(std::vector<Match>& matches, std::vector<Match>& top_matches, int num_peaks) {
	int max_matches = matches.size() > num_peaks ? num_peaks : matches.size();
	while (top_matches.size() < max_matches && matches.size() > 0) {
		auto max_it = matches.begin();
		for (auto it = matches.begin() + 1; it != matches.end(); it++) {
			if (it->val > max_it->val) {
				max_it = it;
			}
		}
		Match n = *max_it;
		matches.erase(max_it);
		// if top_matches doesn't include a AND top_matches doesn't include b
		auto iter_a = std::find_if(top_matches.begin(), top_matches.end(), [&](Match& m){return m.a == n.a;});
		auto iter_b = std::find_if(top_matches.begin(), top_matches.end(), [&](Match& m){return m.b == n.b;});
		if (iter_a == top_matches.end() && iter_b == top_matches.end()) {
			top_matches.push_back(n);
		}
	}
}

void findMatches(FILE* outfile, double threshold, SpectraMGF& spectrum, int num_peaks, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, std::vector<Match>>& peak_matches_shifted, std::unordered_map<std::string, double>& db_precursors, SpectralIndex& spectral_index, bool hybrid_matching) {
	for(auto& [entry_key, precursor]: db_precursors) {
	  std::vector<Match> top_matches;
	  std::vector<Match> top_matches_shifted;

		if (peak_matches[entry_key].size() > 0) {
			getTopMatches(peak_matches[entry_key], top_matches, num_peaks);
		}

		if (!hybrid_matching && peak_matches_shifted[entry_key].size() > 0) {
			getTopMatches(peak_matches_shifted[entry_key], top_matches_shifted, num_peaks);
		}

		// calculate dot product from top matches
		double dot_product = 0;
	  for (int l = 0; l < top_matches.size(); l++) {
	    dot_product += top_matches[l].val;
		//std::cout << "matching peak products:" << top_matches[l].val << std::endl;
	  }
		if (!hybrid_matching) {
			double dot_product_shifted = 0;
			for (int m = 0; m < top_matches_shifted.size(); m++) {
				dot_product_shifted += top_matches_shifted[m].val;
			}
			if (dot_product_shifted > dot_product) {
				dot_product = dot_product_shifted;
				top_matches = top_matches_shifted;
			}
		}
	  if ((dot_product > threshold)) {
	    // Extra file info from meta
	    size_t pos = entry_key.find("#");
	    std::string file_name = spectral_index.getSpectraFile(stoi(entry_key.substr(0, pos)));
	    int scan = stoi(entry_key.substr(pos + 1, std::string::npos));
			double mz_delta = spectrum.mass - db_precursors.at(entry_key);

			//std::cout << entry_key << " precursor: " << precursor <<" product: " << dot_product << " top matches: " << top_matches.size() << std::endl;
			mz_delta = mz_delta < 0 ? -1 * mz_delta : mz_delta;
			fprintf(outfile, "%s\t%d\t%s\t%d\t%lf\t%d\t%lf\n", spectrum.file.c_str(), spectrum.can_num + 1, file_name.c_str(), scan + 1, dot_product, top_matches.size(), mz_delta);
	  }
	}
}

int main(int argc, char** argv){
	std::string query_file_name = "";
  std::string lib_path = "library";
	std::string output_file_name = "matches-all-exact.tsv";
	bool query_is_list = false;
	double threshold=0.2;
  bool error_tolerant = false;
	bool hybrid_matching = false;
	double precursor_tolerance = -1;
  double peak_tolerance = 0.01;
	bool debug = true;
	std::ofstream log_file;
	log_file.open("log_search.txt");
	parse_args(argc, argv, query_file_name, query_is_list, threshold, error_tolerant, hybrid_matching, precursor_tolerance, peak_tolerance, lib_path, output_file_name, debug);

  auto load_start = std::chrono::high_resolution_clock::now();

  //std::cout << "Begin search: " << (error_tolerant ? "error tolerant" : "exact") << ", precursor tolerance: " << precursor_tolerance << ", peak tolerance: " << peak_tolerance << ", threshold: " << threshold << std::endl;

  SpectralIndex spectral_index(lib_path);
  double mass_bin_size = spectral_index.mass_bin_size;
  double pos_bin_size = spectral_index.pos_bin_size;
  //disable peak filtering
  spectral_index.filter_peaks = -1;
  ParseMGF queryMGF (query_file_name, query_is_list, false, spectral_index.filter_peaks, pos_bin_size, mass_bin_size);
	std::vector<std::unordered_map<std::string, std::vector<Match>>> peak_matches{};
	std::vector<std::unordered_map<std::string, std::vector<Match>>> peak_matches_shifted{};
	std::vector<std::unordered_map<std::string, double>> db_precursors{};

	FILE* out_all = fopen(output_file_name.c_str(), "w");
  fprintf(out_all, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Query File", "Query Scan", "DB File", "DB Scan", "Score", "Matched Peaks", "M/Z Delta");


  int parse_status = queryMGF.parseNextFile(log_file);
  if (parse_status < 0) {
		std::cout << "Error parsing " << queryMGF.curr_file_path << std::endl;
    abort();
	}

	auto arbitrary_timestamp = std::chrono::high_resolution_clock::now();
	auto duration_matching = std::chrono::duration_cast<std::chrono::milliseconds>(arbitrary_timestamp - arbitrary_timestamp);
	auto duration_calc = std::chrono::duration_cast<std::chrono::milliseconds>(arbitrary_timestamp - arbitrary_timestamp);
	auto duration_loading = std::chrono::duration_cast<std::chrono::milliseconds>(arbitrary_timestamp - arbitrary_timestamp);
	auto duration_searching = std::chrono::duration_cast<std::chrono::milliseconds>(arbitrary_timestamp - arbitrary_timestamp);

	int num_partitions = spectral_index.num_partitions;
	peak_matches.resize(queryMGF.curr_spectra.size());
	peak_matches_shifted.resize(queryMGF.curr_spectra.size());
	db_precursors.resize(queryMGF.curr_spectra.size());
	for (int part = 0; part < num_partitions; part++) {
		// first load data we will need into memory
	  std::cout << "Loading index..." << std::endl;
	  for (int i = 0; i < queryMGF.curr_spectra.size(); i++) {
	    spectral_index.load(part, queryMGF.curr_spectra[i].mass, queryMGF.curr_spectra[i].peaks, error_tolerant, peak_tolerance);
	  }
	  auto load_stop = std::chrono::high_resolution_clock::now();
	  duration_loading += std::chrono::duration_cast<std::chrono::milliseconds>(load_stop - load_start);
		if (debug) {
			std::cout << "Index loaded. Load time: " << duration_loading.count() << "ms. Searching..." << std::endl;
		} else {
			std::cout << "Index loaded. Searching..." << std::endl;
		}
	  auto search_start = std::chrono::high_resolution_clock::now();

	  duration_matching += std::chrono::duration_cast<std::chrono::milliseconds>(search_start - search_start);
	  duration_calc += std::chrono::duration_cast<std::chrono::milliseconds>(search_start - search_start);

	  // Iterate through spectra in each file
	  for (int i = 0; i < queryMGF.curr_spectra.size(); i++) {
			if (debug) {
				std::cout << "Query spectrum " << i << std::endl;
			}
	    SpectraMGF spectrum = queryMGF.curr_spectra[i];
	    int m = spectrum.mass_rounded;
	    int num_peaks = spectrum.peaks.size();
	    
	    //print the query spectra
		/*
	    std::cout << "BEGIN IONS" << std::endl;
			std::cout << "PEPMASS=" << spectrum.mass << std::endl;
			std::cout << "ROUNDED_MASS=" << spectrum.mass_rounded << std::endl;
			for (int j = 0; j < spectrum.peaks.size(); j++) {
				double pos = spectrum.peaks[j].pos;
				double mag = spectrum.peaks[j].mag;
				std::cout << pos << "\t" << mag << std::endl;
			}
		std::cout << "END IONS" << std::endl;
	    */

	    auto matching_start = std::chrono::high_resolution_clock::now();

	    // Iterate through spectrum peaks
	    for(size_t j = 0; j < num_peaks; j++) {
	      double target_magnitude = spectrum.peaks[j].mag_normalized;
	      int peak_tolerance_units = peak_tolerance / spectral_index.pos_bin_size;
	      for (int t = peak_tolerance_units * -1; t < peak_tolerance_units + 1; t++) {
	        double p = spectrum.peaks[j].pos + t * spectral_index.pos_bin_size;
	        if (error_tolerant) {
	          spectral_index.queryUnshifted(p, j, spectrum.mass, spectrum.peaks[j].pos, target_magnitude, peak_tolerance, peak_matches[i], db_precursors[i]);
	          spectral_index.queryShifted(p, j, spectrum.mass, spectrum.mass - spectrum.peaks[j].pos, target_magnitude, peak_tolerance, hybrid_matching ? peak_matches[i] : peak_matches_shifted[i], db_precursors[i]);
	        } else {
	          spectral_index.queryExact(m, p, spectrum.mass, j, spectrum.peaks[j].pos, target_magnitude, precursor_tolerance, peak_tolerance, peak_matches[i], db_precursors[i]);
	        }
	      }
	    }

	    auto matching_stop = std::chrono::high_resolution_clock::now();
	    duration_matching += std::chrono::duration_cast<std::chrono::milliseconds>(matching_stop - matching_start);
			if (debug) {
				//std::cout << "Finished searching spectrum " << i + 1 << "/" << queryMGF.curr_spectra.size() << "." << std::endl;
			}
	  }
		spectral_index.unload();
		auto search_stop = std::chrono::high_resolution_clock::now();
	  duration_searching += std::chrono::duration_cast<std::chrono::milliseconds>(search_stop - search_start);
	}

	for (int i = 0; i < queryMGF.curr_spectra.size(); i++) {
		auto calc_start = std::chrono::high_resolution_clock::now();
		SpectraMGF spectrum = queryMGF.curr_spectra[i];
		int num_peaks = spectrum.peaks.size();
		findMatches(out_all, threshold, spectrum, num_peaks, peak_matches[i], peak_matches_shifted[i], db_precursors[i], spectral_index, hybrid_matching);
		auto calc_stop = std::chrono::high_resolution_clock::now();
		duration_calc += std::chrono::duration_cast<std::chrono::milliseconds>(calc_stop - calc_start);
	}

	if (debug) {
		std::cout << "Search complete. Search time: " << duration_searching.count() << " ms. " << std::endl;
	} else {
		std::cout << "Search complete." << std::endl;
	}

	if (debug) {
		std::cout << "Match time: " << duration_matching.count() << " ms. " << std::endl;
	  std::cout << "Calc time: " << duration_calc.count() << " ms. " << std::endl;
	}

	fclose(out_all);
	return 0;
}
