#ifndef SPECTRAL_INDEX_INTERFACE_H
#define SPECTRAL_INDEX_INTERFACE_H

#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <sstream>
#include <unordered_map>
#include <cereal/archives/json.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/unordered_map.hpp>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "./mgf/spectra_mgf.hpp"

namespace fs = std::experimental::filesystem;

const std::string null_pos_pointer = "----------------";

struct IndexEntry {
	std::string meta;
	int peak_index;
	double peak_pos;
	double magnitude;
};

struct Match {
	int a;
	int b;
	double val;
};

const unsigned long long max_peaks = 10000000000000000000;

class SpectralIndex{
	public:
		double mass_bin_size;
		double pos_bin_size;
		bool is_new;
		unsigned long long spectra_count = 0;
		unsigned long long peak_count = 0;
		unsigned long long num_partitions = 0;
		int filter_peaks;

		SpectralIndex(std::string library_path, double _mass_bin_size = 0, double _pos_bin_size = 0, double _file_precision = 0.1, int _max_peaks_partition = 100000000, int _filter_peaks = 6);
		void insertSpectrum(int id, std::string filename, int scan, int mass, double mass_raw, int charge, std::vector<Peak>& peaks, std::ofstream& output);
		void commit();
		void load(int part, double precursor, std::vector<Peak>& peaks, bool error_tolerant, double ion_tolerance);
		void unload();
		void queryExact(int precursor, double peak, double query_precursor, int query_peak_index, double query_peak, double query_magnitude, double precursor_tolerance, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list);
		void queryUnshifted(double peak, int query_peak_index, double query_precursor, double query_peak, double query_magnitude, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list);
		void queryShifted(double peak, int query_peak_index, double query_precursor, double query_peak, double query_magnitude, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list);
		std::string getSpectraFile(int ref);

	private:
		double file_precision;
		int max_peaks_partition; // Maximum number of peaks per partition;
		std::string lib_path;
		std::unordered_map<std::string, int> file_ref_map;
		std::vector<std::set<int>> loaded_unshifted;
		std::vector<std::set<int>> loaded_shifted;
		// unshifted_index: partition -> peak file -> peak -> precursor -> data
		std::vector<std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::string>>>> unshifted_index;
		// shifted_index: partition -> difference file (precursor - peak) -> difference -> data
		std::vector<std::unordered_map<int, std::unordered_map<int, std::string>>> shifted_index;

		void findPotentialMatches(std::string entries_str, double query_precursor, int query_peak_index, double query_peak, double query_magnitude, double precursor_tolerance, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list);
		int setSpectraFileRef(std::string filename);
		void updateMeta();
		void addPartition();
};

#endif
