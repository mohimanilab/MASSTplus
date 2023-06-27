#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

std::vector<double> updateNormalizedPeaks(std::vector<int> peak_magnitudes_raw) {
    std::vector<double> peak_magnitudes_normalized;
    double norm = 0;
    for(int i = 0; i < peak_magnitudes_raw.size(); i++){
        norm += (double) peak_magnitudes_raw[i];
    }
    norm = std::sqrt(norm);

    //Normalize vector:
    peak_magnitudes_normalized.resize(peak_magnitudes_raw.size());

    for(int i = 0; i < peak_magnitudes_raw.size(); i++){
        peak_magnitudes_normalized[i] = std::sqrt((double) peak_magnitudes_raw[i]) / norm;
    }

    return peak_magnitudes_normalized;
}

int main(int argc, char** argv){
  std::cout.precision(18);
  std::vector<int> peak_magnitudes_raw;
  peak_magnitudes_raw.push_back(12);
  peak_magnitudes_raw.push_back(37);
  peak_magnitudes_raw.push_back(94);
  peak_magnitudes_raw.push_back(43);
  std::vector<double> peak_magnitudes_normalized = updateNormalizedPeaks(peak_magnitudes_raw);
  for (int i = 0; i < peak_magnitudes_normalized.size(); i++){
    std::cout << peak_magnitudes_normalized[i] << std::endl;
  }
	return 0;
}
