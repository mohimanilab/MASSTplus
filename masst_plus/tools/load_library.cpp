#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../spectral_index.hpp"
#include "../mgf/parse_mgf.hpp"
#include "../mgf/spectra_mgf.hpp"

namespace fs = std::experimental::filesystem;

void print_help_message(){
	printf("Usage: ./load <library_file(list)_name>\n");
	printf("       -h: print this help message\n");
	printf("       --reference-list, -r: indicate library is a file list\n");
	return;
}

void parse_args(int argc, char** argv,
				std::string& library_file_name,
				std::string& library_path,
				bool& library_is_list,
				bool& use_input_stream){

	static struct option long_options[] =
		{
			{"reference-list", no_argument, 0, 'r'},
			{"library", required_argument, 0, 'l'},
			{"use-stream", no_argument, 0, 's'}
		};
	char c;
	int option_index = 0;

	while ((c = getopt_long(argc, argv, "rl:sh", long_options, &option_index)) != -1){
		switch (c){
			case 'r':
				library_is_list = true;
				break;
			case 'l':
				library_path = optarg;
				break;
			case 's':
				use_input_stream = true;
				break;
			case 'h':
				print_help_message();
				exit(0);
			default:
				abort();
		}
	}

	if (argc - optind != 1){
		printf("False number of arguments! \n");
		print_help_message();
		return;
	}
	library_file_name = std::string(argv[optind++]);
	return;
}

int main(int argc, char** argv){
	std::string library_file_name = "";
	bool library_is_list = false;
	bool use_input_stream = false;
	double mass_bin_size = 1.0;
	double pos_bin_size = 0.01;

	double file_precision = 0.1;
	int max_peaks_partition = 100000000;
	const int filter_peaks = -1; // Keep top X peaks in window; -1 = no filtering
	std::string lib_path = "library";
	parse_args(argc, argv, library_file_name, lib_path, library_is_list, use_input_stream);
	SpectralIndex spectral_index (lib_path, mass_bin_size, pos_bin_size, file_precision, max_peaks_partition, filter_peaks);
	std::ofstream log_file;
	log_file.open("log_load.txt");

	std::ofstream spectrum_prints;
	spectrum_prints.open("loaders_prints.mgf");

	auto loading_start = std::chrono::high_resolution_clock::now();
	int spectra_count = 0;
	int total_peaks_count = 0;

	ParseMGF libraryMGF (library_file_name, library_is_list, use_input_stream, filter_peaks, pos_bin_size, mass_bin_size);

	// Iterate through spectra files and load them into index one by one
	if (use_input_stream) {
		std::string stream_file_name = "";
		while(std::getline(std::cin, stream_file_name) && stream_file_name != "END_STREAM") {
			if (stream_file_name == "COMMIT") {
				 log_file << "Saving to file..." << std::endl;
				spectral_index.commit();
				 log_file << "Index saved to file." << std::endl;
			} else {
				log_file << "Parsing " << stream_file_name << std::endl;
				int parse_status = libraryMGF.parseStream(stream_file_name, log_file);
				if (parse_status < 0) {
					 log_file << "Error parsing " << stream_file_name << std::endl;
				} else {
					spectra_count += libraryMGF.curr_spectra.size();
					 log_file << "Loaded " << libraryMGF.curr_spectra.size() << " spectra (" << spectra_count << " total). Inserting..." << std::endl;
					int num_peaks = 0;
					for (int i = 0; i < libraryMGF.curr_spectra.size(); i++) {
						SpectraMGF spectrum = libraryMGF.curr_spectra[i];
						/*
						spectrum_prints << "BEGIN IONS" << std::endl;
						spectrum_prints << "PEPMASS=" << spectrum.mass << std::endl;
						spectrum_prints << "SCAN=" << spectrum.can_num << std::endl;
						spectrum_prints << "ROUNDED_MASS=" << spectrum.mass_rounded << std::endl;
						for (int j = 0; j < spectrum.peaks.size(); j++) {
							double pos = spectrum.peaks[j].pos;
							double mag = spectrum.peaks[j].mag;
							spectrum_prints << pos << "\t" << mag << std::endl;
						}
						spectrum_prints << "END IONS \n" << std::endl;
						*/
						spectral_index.insertSpectrum(spectrum.id, spectrum.file, spectrum.can_num, spectrum.mass_rounded, spectrum.mass, spectrum.charge, spectrum.peaks, log_file);
						num_peaks += spectrum.peaks.size();
					}
					total_peaks_count += num_peaks;
					 log_file << "Inserted " << libraryMGF.curr_spectra.size() << " spectra (" << spectra_count << " total), " << num_peaks << " peaks (" << total_peaks_count << " total)." << std::endl;
				}
			}
		}
	} else {
		while (!libraryMGF.hasFinishedParsing()) {
			int parse_status = libraryMGF.parseNextFile(log_file);
			if (parse_status < 0) {
				 log_file << "Error parsing " << libraryMGF.curr_file_path << std::endl;
			} else {
				spectra_count += libraryMGF.curr_spectra.size();
				 log_file << "Loaded " << libraryMGF.curr_spectra.size() << " spectra (" << spectra_count << " total). Inserting..." << std::endl;
				int num_peaks = 0;
				for (int i = 0; i < libraryMGF.curr_spectra.size(); i++) {
					SpectraMGF spectrum = libraryMGF.curr_spectra[i];
					/*
					spectrum_prints << "BEGIN IONS" << std::endl;
					spectrum_prints << "PEPMASS=" << spectrum.mass << std::endl;
					spectrum_prints << "SCAN=" << spectrum.can_num << std::endl;
					spectrum_prints << "ROUNDED_MASS=" << spectrum.mass_rounded << std::endl;
					for (int j = 0; j < spectrum.peaks.size(); j++) {
						double pos = spectrum.peaks[j].pos;
						double mag = spectrum.peaks[j].mag;
						spectrum_prints << pos << "\t" << mag << std::endl;
					}
					spectrum_prints << "END IONS \n" << std::endl;
					*/
					spectral_index.insertSpectrum(spectrum.id, spectrum.file, spectrum.can_num, spectrum.mass_rounded, spectrum.mass, spectrum.charge, spectrum.peaks, log_file);
					num_peaks += spectrum.peaks.size();
				}
				total_peaks_count += num_peaks;
				 log_file << std::endl << "Inserted " << libraryMGF.curr_spectra.size() << " spectra (" << spectra_count << " total), " << num_peaks << " peaks (" << total_peaks_count << " total)." << std::endl;
				auto current_time = std::chrono::high_resolution_clock::now();
				auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - loading_start);
				 log_file << "Time elapsed: " << time_elapsed.count() << " ms." << std::endl << std::endl;
			}
		}
	}
	spectrum_prints.close();
	 log_file << "Loaded peaks into index. Saving to file..." << std::endl;
	auto saving_start = std::chrono::high_resolution_clock::now();
	spectral_index.commit();
	auto saving_stop = std::chrono::high_resolution_clock::now();
	auto duration_saving = std::chrono::duration_cast<std::chrono::milliseconds>(saving_stop - saving_start);
	 log_file << "Index saved to file." << std::endl;
	if (!use_input_stream) {
		 log_file << "Saving took " << duration_saving.count() << " milliseconds. " << std::endl;
		auto loading_stop = std::chrono::high_resolution_clock::now();
		auto duration_loading = std::chrono::duration_cast<std::chrono::milliseconds>(loading_stop - loading_start);
		 log_file << "Loading took " << duration_loading.count() << " milliseconds. " << std::endl;
		 log_file << "Loaded " << spectra_count << " spectra in the library. " << std::endl;
	} else {
		std::string confirmMessage = "DONE";
		std::string input;
		std::getline(std::cin, input);
		if (input != confirmMessage) {
			log_file << "Error completing script." << std::endl;
		}
	}
	return 0;
}
