#ifndef PARSE_MGF_H9
#define PARSE_MGF_H

#include "spectra_mgf.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <random>

class ParseMGF
{
    public:
      ParseMGF(const std::string& spectra_directory, bool is_list_of_files, bool use_stream, int filter_peaks, double pos_bin_size, double mass_bin_size);
      int parseNextFile(std::ofstream& log_file);
      int parseStream(std::string& spectra_file, std::ofstream& log_file);
      bool hasFinishedParsing();

      double pos_bin_size;
      double mass_bin_size;
      int filter_peaks;
      std::vector <SpectraMGF> curr_spectra;
      std::string curr_file_path;

    private:
      int curr_file_index;
      int curr_spectrum_index;
      std::vector <std::string> spectra_files;

      int parseSpectrum(std::ifstream& file, SpectraMGF& spectrum, const std::string& spectrum_file, int cur_can_num, std::ofstream& log_file);
      int parseSpectrumStream(SpectraMGF& spectrum, const std::string& spectra_file, int cur_can_num, std::ofstream& log_file);

      std::string getFileExtension(const std::string& file_name);
};

#endif
