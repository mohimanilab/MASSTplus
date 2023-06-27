#include "spectra_mgf.hpp"


SpectraMGF::SpectraMGF(){
    file = "";
    can_num = -1;
    return;
}

bool SpectraMGF::operator< (const SpectraMGF& other){
    return (mass_rounded < other.mass_rounded);
}

void SpectraMGF::filterWindowPeaks(int top_window_peaks) {
  const int peak_window_size = 50; // 50 Dalton window in which to select top X peaks
  // Sort peaks
  std::unordered_map<int, std::vector<std::pair<double, double>>> windows{};
  // First, sort into windows
  for (int j = 0; j < peaks.size(); j++) {
    int window = ((int) peaks[j].pos) / peak_window_size;
    if (windows.find(window) == windows.end()) {
      windows.emplace(window, std::vector<std::pair<double, double>>{});
    }
    windows.at(window).push_back(std::make_pair(peaks[j].mag, peaks[j].pos));
  }
  peaks.clear();
  // Iterate through windows and identify top 5 peaks
  for(auto& [_, window_peaks]: windows) {
    // partial sort by peak magnitude
    int sort_peaks = top_window_peaks < window_peaks.size() ? top_window_peaks : window_peaks.size();
    std::partial_sort(window_peaks.begin(), window_peaks.begin()+sort_peaks, window_peaks.end(), std::greater<std::pair<double, double>>());
    for (int k = 0; k < sort_peaks; k++) {
      Peak p;
      p.mag = window_peaks.at(k).first;
      p.pos = window_peaks.at(k).second;
      peaks.push_back(p);
    }
  }
}

// Our normalization
// void SpectraMGF::updateNormalizedPeaks(){
//     float norm = 0;
//     for(int i = 0; i < peak_magnitudes_raw.size(); i++){
//         norm += ((double) peak_magnitudes_raw[i])*((double) peak_magnitudes_raw[i]);
//     }
//     norm = std::sqrt(norm);
//
//     //Normalize vector:
//     peak_magnitudes_normalized.resize(peak_magnitudes_raw.size());
//
//     for(int i = 0; i < peak_magnitudes_raw.size(); i++){
//         peak_magnitudes_normalized[i] = peak_magnitudes_raw[i] / norm;
//     }
//
// }

// MASST normalization
void SpectraMGF::updateNormalizedPeaks() {
    double norm = 0;
    for(int i = 0; i < peaks.size(); i++){
        norm += peaks[i].mag;
    }
    norm = std::sqrt(norm);
    for(int i = 0; i < peaks.size(); i++){
      peaks[i].mag_normalized = std::sqrt(peaks[i].mag) / norm;
    }
}

void SpectraMGF::roundPositions(double interval_size){
    //Get discrete, rounded values for peak positions.
    for(int i = 0; i < peaks.size(); i++) {
        peaks[i].pos_rounded = (int) std::lround((peaks[i].pos) / interval_size);
    }
}

void SpectraMGF::roundMass(double interval_size){
    //Get discrete, rounded value for mass
    mass_rounded = (int) std::lround((mass) / interval_size);
}
