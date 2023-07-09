#include "parse_mgf.hpp"

ParseMGF::ParseMGF(const std::string& spectra_file,
	bool is_list_of_files,
	bool use_stream,
	int filter_peaks,
	double pos_bin_size,
	double mass_bin_size)
: filter_peaks(filter_peaks), pos_bin_size(pos_bin_size), mass_bin_size(mass_bin_size) {
	curr_file_index = 0;

	if (!use_stream) {
		if(!is_list_of_files){
			//Spectra file is a single file
			spectra_files.push_back(spectra_file);
		} else {
			//Spectra file is a list of files, parse all files in the list.
			std::string path;
			std::ifstream file (spectra_file);
			if(file.is_open()){
				while(std::getline(file, path)){
					spectra_files.push_back(path);
				}
			}
		}
	}
}

int ParseMGF::parseNextFile(std::ofstream& log_file){
	curr_spectra.clear();
	// Empty the previous vector from memory
	std::vector<SpectraMGF>().swap(curr_spectra);
	std::string line;
	curr_file_path = spectra_files[curr_file_index];
	// log_file << "Parsing file " << curr_file_index + 1 << "/" << spectra_files.size() << ": " << curr_file_path << "..." << std::endl;
	//std::cout << "Parsing file " << curr_file_index + 1 << "/" << spectra_files.size() << ": " << curr_file_path << "..." << std::endl;

	std::ifstream file (curr_file_path);
	curr_file_index++;

	if(file.is_open()) {
		int cur_can_num = 0;
		while(std::getline(file, line))
		{
			if(line.compare("BEGIN IONS") == 0){
				//Parse each spectrum
				SpectraMGF spectrum;
				int parse_status = parseSpectrum(file, spectrum, curr_file_path, cur_can_num, log_file);
				curr_spectra.push_back(spectrum);
				curr_spectrum_index++;

				if(parse_status < 0){
					file.close();
					return parse_status;
				}
				cur_can_num++;
			}
		}

		file.close();
	} else {
		log_file << "Could not open file: " << curr_file_path << std::endl;
		return -1;
	}
	return 0;
}

bool ParseMGF::hasFinishedParsing(){
	return spectra_files.size() == curr_file_index;
}

int ParseMGF::parseSpectrum(std::ifstream& file, SpectraMGF& spectrum, const std::string& spectra_file, int cur_can_num, std::ofstream& log_file){
	std::string line;

	double mass = -1;
	int charge = 1;

	while(std::getline(file, line) && line.compare("END IONS") != 0) {
		try {
			if (line.substr(0,6).compare("TITLE=") == 0) {
				spectrum.title = line.substr(6);
			} else if (line.substr(0,8).compare("PEPMASS=") == 0) {
				mass = std::stod(line.substr(8));
			} else if (line.substr(0,7).compare("CHARGE=") == 0) {
				try {
					charge = std::stoi(line.substr(7,1));
				} catch(...) {
					// Ignore invalid charge
				}
			} else if (line.substr(0,8).compare("RTINSECONDS=") == 0) {
				continue;
			} else if (line.length() == 0 || !std::isdigit(line[0])) {
				continue;
			} else if (line.find_last_of(" ") != std::string::npos) {
				double pos = std::stod(line.substr(0,line.find_last_of(" ")));
				double mag = std::stod(line.substr(line.find_last_of(" ")+1));
				Peak p;
				p.pos = pos;
				p.mag = mag;
				spectrum.peaks.push_back(p);
			} else {
				continue;
			}
		} catch(...) {
			log_file << "Error parsing " << spectra_file << " " << cur_can_num << std::endl;
		}
	}

	if (!spectrum.title.length()) {
		spectrum.title = spectra_file + "-" + std::to_string(cur_can_num);
	}

	if (mass == -1) {
		std::cout << "Spectrum #" << cur_can_num << " mass missing!" << std::endl;
		return -1;
	}
	spectrum.mass = mass;
	spectrum.charge = charge;
	spectrum.id = curr_spectrum_index;
	curr_spectrum_index++;
	spectrum.file = spectra_file;
	spectrum.can_num = cur_can_num;

	// Filter peaks if necessary
	/*
	if(filter_peaks != -1){
		spectrum.filterWindowPeaks(filter_peaks);
	}
	*/
	// Calculate normalized peak magnitudes and round peak positions, spectrum
	// mass to discrete values
	spectrum.updateNormalizedPeaks();
	spectrum.roundPositions(pos_bin_size);
	spectrum.roundMass(mass_bin_size);

	return 0;
}

std::string ParseMGF::getFileExtension(const std::string& file_name)
{
	if(file_name.find_last_of(".") != std::string::npos)
		return file_name.substr(file_name.find_last_of(".")+1);
	return "";
}

int ParseMGF::parseStream(std::string& spectra_file, std::ofstream& log_file) {
	curr_spectra.clear(); // Empty the previous vector from memory
	log_file << "Parsing mgf from standard input stream: " << spectra_file << std::endl;

	int cur_can_num = 0;
	std::string line;
	std::string end_symbol = "!";
	while(std::getline(std::cin, line) && line != end_symbol) {
		if(line.compare("BEGIN IONS") == 0) {
			//Parse each spectrum
			SpectraMGF spectrum;
			int parse_status = parseSpectrumStream(spectrum, spectra_file, cur_can_num, log_file);

			if(parse_status < 0){
				log_file << "Couldn't parse mgf: " << spectra_file << " " << cur_can_num;
			} else {
				curr_spectra.push_back(spectrum);
			}
			curr_spectrum_index++;
			cur_can_num++;
		}
	}
	return 0;
}

int ParseMGF::parseSpectrumStream(SpectraMGF& spectrum, const std::string& spectra_file, int cur_can_num, std::ofstream& log_file){
	std::string line;
	std::string end_symbol = "!";

	double mass = -1;
	int charge = 1;

	while(std::getline(std::cin, line) && line != end_symbol && line.compare("END IONS") != 0) {
		try {
			if (line.substr(0,6).compare("TITLE=") == 0) {
				spectrum.title = line.substr(6);
			} else if (line.substr(0,8).compare("PEPMASS=") == 0) {
				mass = std::stod(line.substr(8));
			} else if (line.substr(0,7).compare("CHARGE=") == 0) {
				try {
					charge = std::stoi(line.substr(7,1));
				} catch(...) {
					// Ignore invalid charge
				}
			} else if (line.substr(0,8).compare("RTINSECONDS=") == 0) {
				continue;
			} else if (line.length() == 0 || !std::isdigit(line[0])) {
				continue;
			} else if (line.find_last_of(" ") != std::string::npos) {
				double pos = std::stod(line.substr(0,line.find_last_of(" ")));
				double mag = std::stod(line.substr(line.find_last_of(" ")+1));
				Peak p;
				p.pos = pos;
				p.mag = mag;
				spectrum.peaks.push_back(p);
			} else {
				continue;
			}
		} catch (...) {
			log_file << "Error parsing " << spectra_file << " " << cur_can_num << std::endl;
		}
	}

	if (!spectrum.title.length()) {
		spectrum.title = spectra_file + "-" + std::to_string(cur_can_num);
	}

	if (mass == -1) {
		log_file << "Spectrum #" << cur_can_num << " mass missing!" << std::endl;
		return -1;
	}
	spectrum.mass = mass;
	spectrum.charge = charge;
	spectrum.id = curr_spectrum_index;
	curr_spectrum_index++;
	spectrum.file = spectra_file;
	spectrum.can_num = cur_can_num;

	// Filter peaks if necessary
	/*
	if(filter_peaks != -1){
		spectrum.filterWindowPeaks(filter_peaks);
	}
	*/
	// Calculate normalized peak magnitudes and round peak positions, spectrum
	// mass to discrete values
	spectrum.updateNormalizedPeaks();
	spectrum.roundPositions(pos_bin_size);
	spectrum.roundMass(mass_bin_size);

	return 0;
}
