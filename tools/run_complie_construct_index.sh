echo "pll_construct_index"
g++ -g -Wall -Wextra -O3 -Isrc -o pll_construct_index construct_index_main.cc pruned_landmark_labeling.h

echo "h2h_construct_index"
g++ -g -w -Wextra -O3 -Isrc -o h2h_construct_index Tree-Decomposition.index.cpp 
