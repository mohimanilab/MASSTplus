void testMatcher(Matcher& matcher, ParseMGF& parseMGF, int N, std::vector<std::pair<int, int>>& result){
		//Mark beginning of test
		std::cout << std::endl;
		std::cout << "-- Begin of Matcher Test --" << std::endl;

		//Set dataset size
		// std::cout << "============================" << std::endl;
		std::cout << "Testing N = " << N << std::endl;
		parseMGF.setNumSpectra(N);

		std::cout << "Parsed " << parseMGF.spectras.size() << " Spectras" << std::endl;
		//Run matcher
		std::cout << "Start clock" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();

		matcher.match(parseMGF, result);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		std::cout << "Found " << result.size() << " matches in "
				  << duration.count() / 1000.0 << " seconds" << std::endl;

		//Mark ending of test
		std::cout << "-- End of Matcher Test --" << std::endl;
		std::cout << std::endl;
}

void testSearcher(Searcher& searcher, ParseMGF& parseMGF, SpectraMGF& targetMGF, int N, std::vector<int>& result){
		//Mark beginning of test
		// std::cout << std::endl;
		// std::cout << "-- Begin of Searcher Test --" << std::endl;

		//Set dataset size
		// std::cout << "============================" << std::endl;
		// std::cout << "Testing N = " << N << std::endl;
		// parseMGF.setNumSpectra(N);

		// std::cout << "Parsed " << parseMGF.spectras.size() << " Spectras" << std::endl;
		//Run searcher
		// std::cout << "Start clock" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();

		searcher.search(parseMGF, targetMGF, result);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		// std::cout << "Found " << result.size() << " search results in "
				  // << duration.count() / 1000.0 << " seconds" << std::endl;

		//Mark ending of test
		// std::cout << "-- End of Searcher Test --" << std::endl;
		// std::cout << std::endl;
}

void testParams(double res, int filter_peaks, std::string filename, bool is_file, int max_files, int max_spectra, int mode, std::string query_file_name, bool bf){
	auto loading_start = std::chrono::high_resolution_clock::now();
	ParseMGF parseMGF(filename, is_file, filter_peaks, res, 1.0, max_files, max_spectra);
	auto loading_stop = std::chrono::high_resolution_clock::now();
	auto duration_loading = std::chrono::duration_cast<std::chrono::milliseconds>(loading_stop - loading_start);


	//ParseMGF parseMGF("./data", false, filter_peaks, res, 1.0, 25000, 10000000);
	std::cout << "Loading took " << duration_loading.count() << "milliseconds. " << std::endl;
	std::cout << parseMGF.total_num_peaks << std::endl;

	double sparsity = ((double) parseMGF.total_num_peaks) /  (((double)parseMGF.max_pos_index) * (double)(parseMGF.all_spectras.size()) );
	double avg_peaks = ((double) parseMGF.total_num_peaks) /  (double)(parseMGF.all_spectras.size());
	std::cout << "Sparsity Rate: " << sparsity << std::endl;
	std::cout << "Avg peaks per spectra: " << avg_peaks << std::endl;

	std::cout << "max_pos_i: " << parseMGF.max_pos_index << " max_pos: " << parseMGF.max_pos << " min_pos: " << parseMGF.min_pos << std::endl;

	std::cout << "max_mass_i: " << parseMGF.max_mass_index << " max_mass: " << parseMGF.max_mass << " min_mass: " << parseMGF.min_mass << std::endl;

	if (mode == -1 || mode == 0){
		SparseMatcher sparseMatcher(0.8);
		std::vector<std::pair<int, int>> sparseMatcher_result;sparseMatcher_result.clear();

		std::cout << std::endl;
		std::cout << "====================" << std::endl;

		std::cout << "Sparse Matcher: " << std::endl;
		testMatcher(sparseMatcher, parseMGF, -1, sparseMatcher_result); //Match all spectras for sparse matcher

		if (bf){
			BruteforceMatcher bruteforceMatcher(0.8);
			std::vector<std::pair<int, int>> bruteforceMatcher_result;bruteforceMatcher_result.clear();

			std::cout << "--------------------" << std::endl;
			std::cout << "BruteForce Matcher: " << std::endl;
			testMatcher(bruteforceMatcher, parseMGF, -1, bruteforceMatcher_result); //Match all spectras for brute force matcher

			std::sort(sparseMatcher_result.begin(), sparseMatcher_result.end());
			std::sort(bruteforceMatcher_result.begin(), bruteforceMatcher_result.end());

			if (sparseMatcher_result.size() != bruteforceMatcher_result.size()){
				std::cout << "Error: Size doesn't match!" << std::endl;
			}else{
				int unmatched=0;
				for(int i=0; i<bruteforceMatcher_result.size(); i++){
					if (sparseMatcher_result[i] != bruteforceMatcher_result[i]){
						unmatched++;
					}
				}
				if (unmatched){
					std::cout << "Error: " << unmatched << " unmatched results!" << std::endl;
				}else{
				std::cout << "Test passed!" << std::endl;
				}
			}
		}

		std::cout << "====================" << std::endl;
		std::cout << std::endl;
	}

	if (mode == -1 || mode == 1){
		if (query_file_name == ""){
			std::cout << "No query file!" << std::endl;
			return;
		}
		ParseMGF queryMGF(query_file_name, false, filter_peaks, res, 1.0, 1, -1);
		SparseSearcher sparseSearcher(0.8);
		sparseSearcher.construct_searcher(parseMGF);
		std::vector<int> sparseSearcher_result;sparseSearcher_result.clear();

		std::cout << std::endl;
		std::cout << "====================" << std::endl;
		std::cout << "Sparse Searcher: " << std::endl;

		srand(time(0));
		int unmatched = 0;

		auto start = std::chrono::high_resolution_clock::now();
		int cnt=100;
		while(cnt--){
		for(int id = 0; id < 20000; id++){
			sparseSearcher_result.clear();

			//Search a random spectra against others
			testSearcher(sparseSearcher, parseMGF, parseMGF.spectras[id], -1, sparseSearcher_result);
		}
		}
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration_sparse = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		if (bf){
			BruteforceSearcher bruteforceSearcher(0.8);
			std::vector<int> bruteforceSearcher_result;bruteforceSearcher_result.clear();

			start = std::chrono::high_resolution_clock::now();
			for(int id = 0; id < 20000; id++){
				bruteforceSearcher_result.clear();

				// std::cout<<"ID: "<<id<<std::endl;
				// std::cout << "BruteForce Searcher" << std::endl;
				testSearcher(bruteforceSearcher, parseMGF, parseMGF.spectras[id], -1, bruteforceSearcher_result);
				// std::cout<<std::endl;
			}
			stop = std::chrono::high_resolution_clock::now();
			auto duration_bruteforce = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

			if (unmatched){
				std::cout << "Error: " << unmatched << " incorrect searches! " << std::endl;
			}else{
				std::cout << "Searcher all correct! " << std::endl;
			}

			std::cout << "BruteForce searcher took " << duration_bruteforce.count() << " miliseconds. " << std::endl;
		}

		std::cout << "  Sparse   searcher took " << duration_sparse.count() << " miliseconds on 20000 queries for 100 times. " << std::endl;
		std::cout << "====================" << std::endl;
		std::cout << std::endl;
	}



	/*std::vector<int> N = {10000, 100000, 1000000};
	for(int i = 0; i < N.size(); i++){
	std::cout << "====================" << std::endl;
	std::cout << "Sparse Matcher: " << std::endl;
	testMatcher(sparseMatcher, parseMGF, N[i]);

	std::cout << "====================" << std::endl;
	std::cout << "Bruteforce Matcher: " << std::endl;
	testMatcher(bruteforceMatcher, parseMGF, N[i]);
	}*/
}

void print_help_message(){
	std::cout << "Usage: -h: print help message" << std::endl;
	std::cout << "       -i <filename>: specify single input file name" << std::endl;
	std::cout << "       -I <listname>: specify input file list name" << std::endl;
	std::cout << "       -f <max_files>: specify number of files limit, none or -1 for unlimited" << std::endl;
	std::cout << "       -s <max_spectra>: specify number of spectra limit, none or -1 for unlimited" << std::endl;
	std::cout << "       -m <0/1>: 0 to test matcher, 1 to test searcher" << std::endl;
	std::cout << "       -q <filename>: specify query file in searcher mode" << std::endl;
	std::cout << "       -b: switch on tests of bruteforce version" << std::endl;
	std::cout << "       -e: Run an example" << std::endl;
	std::cout << std::endl;
	return;
}

void parse_args(int argc, char** argv, std::string& library_file_name, bool& is_file, int& max_files, int& max_spectra, int& mode, std::string& query_file_name, bool& bf){
	bool has_input_library = false;
	bool has_arguments = false;
	char c = '\0';
	while((c = getopt(argc, argv, "hi:I:f:s:m:q:be")) != -1){
		switch(c){
			case 'h':
				// std::cout <<"Got h" << std::endl;
				print_help_message();
				exit(0);
			case 'i':
				// std::cout <<"Got i" << std::endl;
				if (has_input_library){
					std::cerr << "Error: Please specify either input file or input directory! " << std::endl;
					exit(0);
				}
				has_input_library = true;
				library_file_name = (std::string) optarg;
				is_file = true;
				break;
			case 'I':
				// std::cout <<"Got I" << std::endl;
				if (has_input_library){
					std::cout << "Error: Please specify either input file or input directory! " << std::endl;
					exit(0);
				}
				has_input_library = true;
				library_file_name = (std::string) optarg;
				is_file = false;
				break;
			case 'f':
				// std::cout <<"Got f" << std::endl;
				max_files = atoi(optarg);
				break;
			case 's':
				// std::cout <<"Got s" << std::endl;
				max_spectra = atoi(optarg);
				break;
			case 'm':
				mode = atoi(optarg);
				break;
			case 'q':
				query_file_name = (std::string) optarg;
			case 'e':
				// std::cout <<"Got e" << std::endl;
				if (!has_arguments){
					std::cout << "Running examples..." << std::endl;
					max_files = 100;
					max_spectra = -1;
					return;
				}
				break;
			case 'b':
				bf = true;
			default:
				// std::cout <<"?" << std::endl;
				std::cerr << "Error: No such argument! " << std::endl;
				print_help_message();
				exit(0);
		}
		has_arguments = true;
	}
	if (!has_arguments){
		print_help_message();
		exit(0);
	}
	return;
}

int main(int argc, char** argv) {
	// freopen("output.txt", "w", stdout);
	int max_files=-1, max_spectra=-1;
	std::string library_file_name = "/projects/mohimanilab/ben/sparsematch/data/";
	std::string query_file_name = "";
	int mode = -1;
	bool is_file = false;
	bool bf = false;
	parse_args(argc, argv, library_file_name, is_file, max_files, max_spectra, mode, query_file_name, bf);
	testParams(0.01, 5, library_file_name, is_file, max_files, max_spectra, mode, query_file_name, bf);
	/*std::vector<int> filter = {5, 10, 20, 50, 100};
	for(int i = 0; i < filter.size(); i++){
	std::cout << "======================================" << std::endl;
	std::cout << "Filter test: " << filter[i] << std::endl;
	std::cout << "======================================" << std::endl;
		testParams(0.01, filter[i]);
	}*/

	/*std::vector<double> res = {0.01, 0.02, 0.05, 0.1};
	for(int i = 0; i < res.size(); i++){
		std::cout << "======================================" << std::endl;
	std::cout << "Resolution test: " << res[i] << std::endl;
	std::cout << "======================================" << std::endl;
	testParams(res[i], -1);
	}*/
	// fclose(stdout);
}

