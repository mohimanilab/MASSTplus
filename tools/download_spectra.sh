echo 'Downloading all spectra files...'
wget -c -np -nc -P data/ ftp://massive.ucsd.edu/MSV000080673/ccms_peak/2017.AmericanGut3K.mzXMLfiles/Samples/*
echo 'Finished downloading spectra files'
