#include "spectral_index.hpp"

SpectralIndex::SpectralIndex(std::string library_path, double _mass_bin_size, double _pos_bin_size, double _file_precision, int _max_peaks_partition, int _filter_peaks) {
	lib_path = library_path;
	if (!_mass_bin_size || !_pos_bin_size || fs::exists(lib_path + "/config.txt")) {
		is_new = false;
		// Use existing library
		FILE* config = fopen((lib_path + "/config.txt").c_str(), "r");
		if (!config){
			std::cout << "Library corrupted!" << std::endl;
			return;
		}
		fscanf(config, "%lf %lf %lf %d %d", &mass_bin_size, &pos_bin_size, &file_precision, &max_peaks_partition, &filter_peaks);
		fclose(config);
		// Load metadata
		FILE* meta = fopen((lib_path + "/meta.txt").c_str(), "r");
		if (!meta) {
			std::cout << "Unable to open meta file" << std::endl;
			return;
		}
		fscanf(meta, "%d %d %d", &spectra_count, &peak_count, &num_partitions);
		fclose(meta);
		loaded_unshifted.resize(num_partitions);
		loaded_shifted.resize(num_partitions);
		unshifted_index.resize(num_partitions);
		shifted_index.resize(num_partitions);
		// Load spectra file refs
		std::ifstream ifile;
		ifile.open((lib_path + "/spectra_ref.txt").c_str());
		cereal::JSONInputArchive iarchive(ifile);
		iarchive(file_ref_map);
	} else {
		is_new = true;
		// Initialize library
		mass_bin_size = _mass_bin_size;
		pos_bin_size = _pos_bin_size;
		file_precision = _file_precision;
		max_peaks_partition = _max_peaks_partition;
		filter_peaks = _filter_peaks;
		// Create library directory
		if (!fs::exists(lib_path)) {
			int check = mkdir(lib_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if (check == -1) {
				std::cout << "Failed to create library directory! " << std::endl;
				return;
			}
		}
		// Create config file
		FILE* config = fopen((lib_path + "/config.txt").c_str(), "w");
		if (!config){
			std::cout << "Failed to create configuration file! " << std::endl;
		}
		fprintf(config, "%lf %lf %lf %d %d", mass_bin_size, pos_bin_size, file_precision, max_peaks_partition, filter_peaks);
		fclose(config);
		// Create first partition
		mkdir((lib_path + "/partitions").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		addPartition();
		// Create metadata file
		updateMeta();
	}
	return;
}

void SpectralIndex::insertSpectrum(int id, std::string filename, int scan, int mass, double mass_raw, int charge, std::vector<Peak>& peaks, std::ofstream& output) {
	if (!is_new) {
		// load any existing entries so we don't overwrite them
		int current_part = num_partitions - 1;
		load(current_part, mass_raw, peaks, false, 0);
	}
	int fileref = file_ref_map.count(filename) ? file_ref_map[filename] : setSpectraFileRef(filename);
	spectra_count++;
	if (peak_count > num_partitions * max_peaks_partition) {
		// Time to move to a new partition
		output << "Saving to file..." << std::endl;
		commit();
		output << "Index saved to file." << std::endl;
		addPartition();
		unshifted_index.at(num_partitions - 2).clear();
		shifted_index.at(num_partitions - 2).clear();
		output << "Moved to new partition. #partitions: " << num_partitions << std::endl;
		if (spectra_count > max_peaks) {
			// This is a safe point at which to terminate
			output << "Maximum number of peaks reached." << num_partitions << std::endl;
			abort();
		}
	}
	for (int i = 0; i < peaks.size(); i++) {
		peak_count++;
		std::string entry_unshifted = std::to_string(fileref) + "#" + std::to_string(scan) + "," + std::to_string(mass_raw) + "," + std::to_string(i) + "," + std::to_string(peaks[i].pos) + "," + std::to_string(peaks[i].mag_normalized);
		std::string entry_shifted = std::to_string(fileref) + "#" + std::to_string(scan) + "," + std::to_string(mass_raw) + "," + std::to_string(i) + "," + std::to_string(mass_raw * charge - peaks[i].pos) + "," + std::to_string(peaks[i].mag_normalized);
		// update unshifted index
		int peak_file = std::round(peaks[i].pos_rounded * file_precision);
		if (unshifted_index.at(num_partitions - 1).find(peak_file) == unshifted_index.at(num_partitions - 1).end()) {
			unshifted_index.at(num_partitions - 1).emplace(peak_file, std::unordered_map<int, std::unordered_map<int, std::string>>{});
		}
		if (unshifted_index.at(num_partitions - 1)[peak_file].find(peaks[i].pos_rounded) == unshifted_index.at(num_partitions - 1)[peak_file].end()) {
			unshifted_index.at(num_partitions - 1)[peak_file].emplace(peaks[i].pos_rounded, std::unordered_map<int, std::string>{});
			unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded].emplace(mass, entry_unshifted);
		} else {
			if (unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded].find(mass) == unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded].end()) {
				unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded].emplace(mass, entry_unshifted);
			} else {
				unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded][mass] = unshifted_index.at(num_partitions - 1)[peak_file][peaks[i].pos_rounded][mass] + "|" + entry_unshifted;
			}
		}
		// update shifted index
		int diff = std::round((mass_raw * charge - peaks[i].pos) / pos_bin_size);
		int diff_file = std::round(diff * file_precision);
		if (diff_file < 0) {
			continue;
		}
		if (shifted_index.at(num_partitions - 1).find(diff_file) == shifted_index.at(num_partitions - 1).end()) {
			shifted_index.at(num_partitions - 1).emplace(diff_file, std::unordered_map<int, std::string>{});
		}
		if (shifted_index.at(num_partitions - 1)[diff_file].find(diff) == shifted_index.at(num_partitions - 1)[diff_file].end()) {
			shifted_index.at(num_partitions - 1)[diff_file].emplace(diff, entry_shifted);
		} else {
			shifted_index.at(num_partitions - 1)[diff_file][diff] = shifted_index.at(num_partitions - 1)[diff_file][diff] + "|" + entry_shifted;
		}
	}
}

// Writes data from queue to database file and header map back to file;
// until this is called, inserted data has not actually been persisted to db.
// Should call after doing a batch of inserts.
void SpectralIndex::commit() {
	// Write unshifted index files
	mkdir((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/unshifted").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	for (auto& [peak_file, unshifted_data]: unshifted_index.at(num_partitions - 1)) {
		std::ofstream unshifted_ofile;
		unshifted_ofile.open((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/unshifted/" + std::to_string(peak_file) + ".txt").c_str());
		cereal::JSONOutputArchive unshifted_archive(unshifted_ofile);
		unshifted_archive(unshifted_data);
	}
	mkdir((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/shifted").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	for (auto& [difference_file, shifted_data]: shifted_index.at(num_partitions - 1)) {
		std::ofstream shifted_ofile;
		shifted_ofile.open((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/shifted/" + std::to_string(difference_file) + ".txt").c_str());
		cereal::JSONOutputArchive shifted_archive(shifted_ofile);
		shifted_archive(shifted_data);
	}
	updateMeta();
}

void SpectralIndex::load(int part, double precursor, std::vector<Peak>& peaks, bool error_tolerant, double ion_tolerance) {
	for(int i = 0; i < peaks.size(); i++) {
		// Unshifted index
		int min_peak = std::round(((peaks[i].pos - ion_tolerance) / pos_bin_size) * file_precision);
		int max_peak = std::round(((peaks[i].pos + ion_tolerance) / pos_bin_size) * file_precision);
		for (int peak_file = min_peak; peak_file <= max_peak; peak_file++) {
			if (loaded_unshifted.at(part).find(peak_file) == loaded_unshifted.at(part).end()) {
				loaded_unshifted.at(part).insert(peak_file);
				std::unordered_map<int, std::unordered_map<int, std::string>> unshifted_part;
				std::ifstream unshifted_file;
				unshifted_file.open((lib_path + "/partitions/partition-" + std::to_string(part) + "/unshifted/" + std::to_string(peak_file) + ".txt").c_str());
				if (unshifted_file.is_open()) {
					cereal::JSONInputArchive unshifted_iarchive(unshifted_file);
					unshifted_iarchive(unshifted_part);
					unshifted_index.at(part)[peak_file].insert(unshifted_part.begin(), unshifted_part.end());
				}
			}
		}
		// Shifted index
		if (error_tolerant) {
			int min_peak_bin = std::round((precursor - peaks[i].pos - ion_tolerance) / pos_bin_size);
			int max_peak_bin = std::round((precursor - peaks[i].pos + ion_tolerance) / pos_bin_size);
			int min_diff = std::round(min_peak_bin * file_precision);
			int max_diff = std::round(max_peak_bin * file_precision);
			for (int diff = min_diff; diff <= max_diff; diff++) {
				if (loaded_shifted.at(part).find(diff) == loaded_shifted.at(part).end()) {
					loaded_shifted.at(part).insert(diff);
					std::unordered_map<int, std::string> shifted_part;
					std::ifstream shifted_file;
					shifted_file.open((lib_path + "/partitions/partition-" + std::to_string(part) + "/shifted/" + std::to_string(diff) + ".txt").c_str());
					if (shifted_file.is_open()) {
						cereal::JSONInputArchive shifted_iarchive(shifted_file);
						shifted_iarchive(shifted_part);
						shifted_index.at(part)[diff].insert(shifted_part.begin(), shifted_part.end());
					}
				}
			}
		}
	}
}

// Free up memory
void SpectralIndex::unload() {
	unshifted_index.clear();
	shifted_index.clear();
	loaded_unshifted.clear();
	loaded_shifted.clear();
	loaded_unshifted.resize(num_partitions);
	loaded_shifted.resize(num_partitions);
	unshifted_index.resize(num_partitions);
	shifted_index.resize(num_partitions);
}

void SpectralIndex::queryExact(int precursor, double peak, double query_precursor, int query_peak_index, double query_peak, double query_magnitude, double precursor_tolerance, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list) {
	int peak_bin = std::round(peak / pos_bin_size);
	int peak_file = std::round(peak_bin * file_precision);
	for (int part = 0; part < num_partitions; part++) {
		if (unshifted_index.at(part).find(peak_file) != unshifted_index.at(part).end()) {
			if (unshifted_index.at(part)[peak_file].find(peak_bin) != unshifted_index.at(part)[peak_file].end() && unshifted_index.at(part)[peak_file][peak_bin].find(precursor) != unshifted_index.at(part)[peak_file][peak_bin].end()) {
				findPotentialMatches(unshifted_index.at(part)[peak_file][peak_bin][precursor], query_precursor, query_peak_index, query_peak, query_magnitude, precursor_tolerance, peak_tolerance, peak_matches, precursors_list);
			}
		}
	}
}

void SpectralIndex::queryUnshifted(double peak, int query_peak_index, double query_precursor, double query_peak, double query_magnitude, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list) {
	int peak_bin = std::round(peak / pos_bin_size);
	int peak_file = std::round(peak_bin * file_precision);
	for (int part = 0; part < num_partitions; part++) {
		if (unshifted_index.at(part).find(peak_file) != unshifted_index.at(part).end() && unshifted_index.at(part)[peak_file].find(peak_bin) != unshifted_index.at(part)[peak_file].end()) {
			for (auto& [precursor, entries_str] : unshifted_index.at(part)[peak_file][peak_bin]) {
				findPotentialMatches(entries_str, query_precursor, query_peak_index, query_peak, query_magnitude, -1, peak_tolerance, peak_matches, precursors_list);
			}
		}
	}
}

void SpectralIndex::queryShifted(double peak, int query_peak_index, double query_precursor, double query_peak, double query_magnitude, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list) {
	int diff = std::round((query_precursor - peak) / pos_bin_size);
	int diff_file = std::round(diff * file_precision);
	for (int part = 0; part < num_partitions; part++) {
		if (shifted_index.at(part).find(diff_file) != shifted_index.at(part).end() && shifted_index.at(part)[diff_file].find(diff) != shifted_index.at(part)[diff_file].end()) {
			findPotentialMatches(shifted_index.at(part)[diff_file][diff], query_precursor, query_peak_index, query_peak, query_magnitude, -1, peak_tolerance, peak_matches, precursors_list);
		}
	}
}

void SpectralIndex::findPotentialMatches(std::string entries_str, double query_precursor, int query_peak_index, double query_peak, double query_magnitude, double precursor_tolerance, double peak_tolerance, std::unordered_map<std::string, std::vector<Match>>& peak_matches, std::unordered_map<std::string, double>& precursors_list) {
	// First split string into multiple entries
	size_t pos = 0;
	bool done = false;
	double db_precursor;
	while (!done) {
		std::string entry_str;
		pos = entries_str.find("|");
		if (pos == std::string::npos) {
			entry_str = entries_str;
			done = true;
		} else {
			entry_str = entries_str.substr(0, pos);
			entries_str.erase(0, pos + 1);
		}

		// Extract data from entry string
		IndexEntry entry;
		size_t pos_e = 0;
		for (int i = 0; i < 4; i++) {
			pos_e = entry_str.find(",");
			switch(i) {
				case 0:
					entry.meta = entry_str.substr(0, pos_e);
					break;
				case 1:
					db_precursor = std::stod(entry_str.substr(0, pos_e));
					break;
				case 2:
					entry.peak_index = std::stoi(entry_str.substr(0, pos_e));
					break;
				case 3:
					entry.peak_pos = std::stod(entry_str.substr(0, pos_e));
					break;
			}
			entry_str.erase(0, pos_e + 1);
		}
		entry.magnitude = std::stod(entry_str);

		if ((precursor_tolerance < 0 || std::abs(db_precursor - query_precursor) <= precursor_tolerance) && std::abs(query_peak - entry.peak_pos) <= peak_tolerance) {
			// Update peak_matches
			if (peak_matches.find(entry.meta) == peak_matches.end()) {
				peak_matches.emplace(entry.meta, std::vector<Match>{});
			}
			if (precursors_list.find(entry.meta) == precursors_list.end()) {
				precursors_list.emplace(entry.meta, db_precursor);
			}
			Match match;
			match.a = query_peak_index;
			match.b = entry.peak_index;
			match.val = entry.magnitude * query_magnitude;
			peak_matches[entry.meta].push_back(match);
			std::cout << "peak found: " << entry.meta << " " << entry.peak_pos << " " << entry.magnitude << std::endl;
		}
	}
}

// Converts a reference file path to an integer that represents it in the index
int SpectralIndex::setSpectraFileRef(std::string filename) {
	std::hash<std::string> hasher;
	int hashed = hasher(filename);
	file_ref_map[filename] = std::abs(hashed);
	std::ofstream ofile;
	ofile.open((lib_path + "/spectra_ref.txt").c_str());
	cereal::JSONOutputArchive file_ref_archive(ofile);
	file_ref_archive(file_ref_map);
	return file_ref_map[filename];
}

// Returns filename corresponding to a given index id and file reference
std::string SpectralIndex::getSpectraFile(int ref) {
	std::string filename = "";
	for (auto& [key, val]: file_ref_map) {
		if (val == ref) filename = key;
	}
	return filename.length() ? filename : "file-not-found";
}

void SpectralIndex::updateMeta() {
	FILE* meta = fopen((lib_path + "/meta.txt").c_str(), "w");
	if (!meta){
		std::cout << "Failed to write to meta file! " << std::endl;
	}
	fprintf(meta, "%llu %llu %llu", spectra_count, peak_count, num_partitions);
	fclose(meta);
}

void SpectralIndex::addPartition() {
	num_partitions++;
	mkdir((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/shifted").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir((lib_path + "/partitions/partition-" + std::to_string(num_partitions - 1) + "/unshifted").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	loaded_unshifted.resize(num_partitions);
	loaded_shifted.resize(num_partitions);
	unshifted_index.resize(num_partitions);
	shifted_index.resize(num_partitions);
	is_new = true;
	updateMeta();
}
