#ifndef SPECTRA_MGF_H
#define SPECTRA_MGF_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <climits>
#include <cmath>

struct Peak {
	double pos; //position
  int pos_rounded;
  double mag; // magnitude
  double mag_normalized;
};

class SpectraMGF
{
    public:
        std::string file;
        int can_num;
        std::string title;
        int id;

        std::vector<Peak> peaks;

        double mass;
        int mass_rounded;
        int charge;


        SpectraMGF();
        void filterWindowPeaks(int top_window_peaks);
        void updateNormalizedPeaks();

        void roundPositions(double interval_size);
        void roundMass(double interval_size);

        bool operator< (const SpectraMGF& other);
 };


#endif
