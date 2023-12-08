#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
using namespace std;
int readGraphHeader(char* file){
    string ifstr = file;
    ifstream ifs;
    ifs.open(ifstr);
    int n;
    ifs >> n;
    cout<<"node# "<<n<<endl;
    ifs.close();
    return n;
}
string remove_gr(char* file){
    string ingr = file;
    size_t pos = ingr.find(".gr");

    // Check if the suffix was found
    if (pos != std::string::npos) {
        // Remove the suffix
        ingr = ingr.substr(0, pos);
    }

    std::cout << ingr << std::endl; // Output: NY
    return ingr;
}
// in param
// ./rdgen.out NY.gr 11 50 
//                          expected output
//                          NY_11.Q
int main(int argc, char** argv) {
    //in params and outfile
    for(int i=1; i<4; i++){
        cout<<argv[i]<<" ";
    }
    cout<<endl;
    int n=0;
    n = readGraphHeader(argv[1]);

    string cityName = remove_gr(argv[1]);
    string qsize = argv[2];
    int querySize = stoi(argv[2]);
    cout<<"qSize: "<<querySize<<endl;
    int qSetNum = stoi(argv[3]);
    cout<<"#Q set: "<<qSetNum<<endl;

    string outQName = cityName+"_"+qsize+".Q";
    cout<<"out: "<<outQName<<endl;

    //gen rand
    std::random_device rd;
    std::mt19937 gen(rd());

    vector<int> shuffle_vec;
    for(int i=1; i<=n; i++){
        shuffle_vec.push_back(i);
    }
    shuffle(shuffle_vec.begin(), shuffle_vec.end(), gen);

    unordered_set<int> usi;
    std::uniform_int_distribution<int> dist(1, n);

    while(usi.size() < qSetNum){
        int randomNumber = dist(gen);
        usi.emplace(randomNumber);
    }
    //cargo Q file
    vector<vector<int>> outQvec;
    for(auto ele:usi){
        vector<int> tmp;
        int startPos = ele;
        //11
        for(int i=0; i<querySize; i++){
            int pos = (startPos+i)%(n+1);
            tmp.push_back(shuffle_vec[pos]);
        }
        outQvec.push_back(tmp);
    }
    //out to qfile
    ofstream ofs;
    ofs.open(outQName);
    for(auto line:outQvec){
        for(auto ele:line){
            ofs<<ele<<" ";
        }
        ofs<<endl;
    }
    ofs.close();
    return 0; 
}
