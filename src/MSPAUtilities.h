#pragma once
#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <queue>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <xmmintrin.h>
#include <random>
#include "H2H.h"
#include "pruned_landmark_labeling.h"
#include "appx.h"
using namespace std;

#define LLINF 9223372036854775807
#define INF 2147483647 
#define TRIQUATINF 1431655765
#define MYMAX_INT 2147483647
#define MYMAX_INT_HALF 1073741823
#define NON -1
// #define RECOVER_MSPATREE

char* outfile;
FILE* outOpenFile;

char* incandfile;
FILE* incandOpenFile;

//res2File
vector<vector<long long>> res50Dist(6);
vector<vector<double>> res50Time(6);


//RECOVER_MSPATREE
struct prevTriState{
    int root;
    int bgt1;
    int msk1;
    int bgt2;
    int msk2;
    prevTriState(){}
    prevTriState(int r, int b1, int m1, int b2, int m2):\
    root(r), bgt1(b1), msk1(m1), bgt2(b2), msk2(m2){}
};
struct prevState{
    int root;
    int bgt1;
    int msk1;
    int apflag1;
    int bgt2;
    int msk2;
    int apflag2;
    prevState(){}
    prevState(int r, int b1, int m1, int f1, int b2, int m2, int f2):\
    root(r), bgt1(b1), msk1(m1), apflag1(f1), bgt2(b2), msk2(m2), apflag2(f2){}
};


template<typename T>
struct ComparePairs{
    bool operator()(pair<int, T>& pa, pair<int, T>& pb){
        return pa.first<pb.first;
    }
};

// apx
struct Spliter{
    vector<int> assignSrcSet;
    int subtreeRoot;
    int prevRoot;
    int subtreeWeight;
    Spliter(){}
    Spliter(vector<int>& inVec, int subR, int prevR, int subWei):\
    assignSrcSet(inVec), subtreeRoot(subR), prevRoot(prevR), subtreeWeight(subWei){}
};
//return "1" bits in a num
//(0011)==3, bitCnt(3) return 2
int bitCnt(int mask){
    bitset<64> b(mask);
    return b.count();
}
int bitDist( vector<int>&terminals, int mask, int nodeid){
    int dis = 0;
    bitset<64> bs(mask);
    for(int i = 0; i<(int)bs.size(); i++){
        if(1 == bs[i]){
            dis += distanceQuery(nodeid, terminals[i]);
        }
    }
    return dis;
}
int char2int(char* pc){
    int i=0;
    while(pc[i]!='\0'){
        i++;
    }
    int res=0;
    for(int cnt=i-1; cnt>=0; cnt--){
        res+=((pc[cnt]-'0')*pow(10, (i-cnt-1)));
    }
    return res;
}
//retrieve assembly points for basic dp
void dfsRTree_BasicDP(prevTriState*** triTree, int root, int rcnt, int mask){
    if(0 == root) return;
    if(0 == rcnt && 0 == mask) return;
    if(rcnt >=1 && mask >0){
        cout<<"assembly point: "<<root<<" budget "<<rcnt<<" mask: "<<mask<<endl;
    }
    prevTriState pts = triTree[mask][root][rcnt];
    int r = pts.root;
    int bgt1 = pts.bgt1;
    int mask1 = pts.msk1;
    int bgt2 = pts.bgt2;
    int mask2 = pts.msk2;
    dfsRTree_BasicDP(triTree, r, bgt1, mask1);
    dfsRTree_BasicDP(triTree, r, bgt2, mask2);
}
//retrieve assembly points for adv dp
void dfsRetrieveTree(prevState**** rTree, int root, int rcnt, int mask, int apflag){
    if(0 == root) return;
    if(0 == rcnt && 0 == mask && 0 == apflag) return; 
    if(rcnt >= 1 && 1 == apflag && mask > 0){
        cout<<"assembly point: "<<root<<" budget: "<<rcnt<<" mask: "<<mask<<endl;
    }
    prevState pst = rTree[mask][root][rcnt][apflag];
    int r = pst.root;
    int bgt1 = pst.bgt1;
    int mask1 = pst.msk1;
    int apflag1 = pst.apflag1;
    int bgt2 = pst.bgt2;
    int mask2 = pst.msk2;
    int apflag2 = pst.apflag2;

    dfsRetrieveTree(rTree, r, bgt1, mask1, apflag1);
    dfsRetrieveTree(rTree, r, bgt2, mask2, apflag2);
}

void loadTerminals(vector<vector<int>>& terminals, char* queryFile, int srcNsink_num){
    FILE* pQf = fopen(queryFile, "r");
    vector<int> Q;
    int tmpNodeID = 0;
    while(EOF != fscanf(pQf, "%d", &tmpNodeID)){
        Q.push_back(tmpNodeID);
    }
    vector<int>  tmp;
    for(int j=0; j<(int)Q.size(); j++){
        tmp.push_back(Q[j]);
        if((int)tmp.size() == srcNsink_num){
            terminals.push_back(tmp);
            tmp.clear();
            continue;
        }
    }
    fclose(pQf);
}
void printRes2File(){
    for(int i=1; i<(int)res50Dist.size(); i++){
        if(res50Dist[i].size() == 0) {
            // cout<<i<<"th res is void"<<endl;
            continue;
        }
        // cout<<i<<"th size "<<res50Dist[i].size()<<endl; 
        for(int j=0; j<(int)res50Dist[i].size(); j++){
            fprintf(outOpenFile, "%ld", res50Dist[i][j]); 
            putc('\t', outOpenFile);
            fprintf(outOpenFile, "%lf", res50Time[i][j]); 
            putc('\n', outOpenFile);
        }
        double tspan = 0.0;

        // cout<<i<<"th avg dist and time"<<endl; 
        int invalid = 0;
        for(int k=0; k<(int)res50Dist[i].size(); k++){
            if(res50Time[i][k] > MYMAX_INT_HALF){
                invalid++;
                continue;
            }
            tspan+=res50Time[i][k];
        }
        int valid = res50Dist[i].size() - invalid;
        if(0 == valid){
            tspan = INF;
        }else{
            tspan/=valid;
        }
        fprintf(outOpenFile, "%d%s", i, "th Avg Time " ); 
        putc('\t', outOpenFile);
        fprintf(outOpenFile, "%lf", tspan); 
        putc('\n', outOpenFile);
    }
}
void printCand2File(string icf, vector<vector<int>>& cand){
    fstream fs;
    fs.open(icf, ios::out);
    for(auto line : cand){
        for(int i=0; i<line.size(); i++){
            if(i+1 == line.size()){
                fs<<line[i]<<'\n';
            }
            fs<<line[i]<<" ";
        }
    }
    fs.close();
}
void fillG(vector<int>& src, int cand, vector<vector<pair<int, int>>>& newG){
    // src
    for (int i = 0; i < src.size(); i++)
    {
        for (int j = 0; j < src.size(); j++)
        {
            newG[i][j] = make_pair(j, distanceQuery(src[i], src[j]));
        }
        int newCur = src.size();
        newG[i][newCur] = make_pair(newCur, distanceQuery(src[i], cand));
    }
    // aps
    int apsCur = src.size();
    for (int j = 0; j < src.size(); j++)
    {
        newG[apsCur][j] = make_pair(j, distanceQuery(cand, src[j]));
    }
    // int newCur = src.size();
    newG[apsCur][apsCur] = make_pair(apsCur, distanceQuery(cand, cand));
    // for(int i=0; i<newG.size(); i++){
    //     cout<<i<<endl;
    //     for(int j=0; j<newG[i].size(); j++){
    //         cout<<newG[i][j].second<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
}
int lowerBound(vector<int>& apset, vector<int>& terminal){
    int lb=0;
    vector<int> apxSteinerTerminals(apset.begin(), apset.end());
    for(auto ele:terminal){
        apxSteinerTerminals.push_back(ele);
    }
    vector<vector<pair<int, int>>> newG;
    fillG(newG, apxSteinerTerminals);
    lb = prim(newG, apxSteinerTerminals.size());
    lb = lb/2;
    return lb;
}

std::chrono::time_point<std::chrono::high_resolution_clock> timerSt(){
   return std::chrono::high_resolution_clock::now();
}
std::chrono::time_point<std::chrono::high_resolution_clock> timerNd(){
   return std::chrono::high_resolution_clock::now();
}
double dur(chrono::time_point<std::chrono::high_resolution_clock> nd, chrono::time_point<std::chrono::high_resolution_clock> st){
    std::chrono::microseconds tot_time = \
    std::chrono::duration_cast<std::chrono::microseconds>(nd - st);
    return (double)(1.0*tot_time.count()/1000000);
}
//We set each time out as half an hour such that 50 queries can finish around a day.
int timeOutCheck(std::chrono::time_point<std::chrono::high_resolution_clock> algo_st, double timeOutBound=1800){
    int flag = 0;
    auto checkTimer = timerNd();
    double timeout = dur(checkTimer, algo_st);
    if (timeout >= timeOutBound){
        flag = 1;
    }
    return flag;
}
// //algo2 with idx, no path retrieval.
// void algo2_ComputePathForSource_with_idx(vector<int>& terminal, vector<int>& cand, \
// unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
// int kth, int& s2c){
//     int sink = terminal.back();
//     vector<int> src(terminal.begin(), terminal.end()-1);
//     //will be delivered to algo3, and will negelect the unchosen cand
//     for(int i=0; i<kth; i++){
//         umpii.emplace(make_pair(cand[i], 0));
//     }
//     //allocate to the nearest cand
//     for(int t=0; t<src.size(); t++){
//         int tmpSrc = src[t];
//         int src2sink = distanceQuery(t, sink);
//         int sc2cd = INF;
//         int chosenCand = 0;
//         //choose cand limited by budget
//         for(int i=0; i<kth; i++){
//             int tmpSrc2cand = distanceQuery(tmpSrc, cand[i]);
//             if(sc2cd > tmpSrc2cand){
//                 sc2cd = tmpSrc2cand; 
//                 chosenCand = i;
//             }
//         }
//         if(sc2cd>src2sink){
//             s2c+=src2sink;
//         // if not assigned to sink, then mark cand chosen by one src.
//         }else{
//             s2c+=sc2cd;
//             umpii[cand[chosenCand]] += 1;
//         }
//     }
//     // cout<<"s2c:"<<s2c<<endl;
// }
// //algo3, no path retrieval.
// void algo3_computePathAssembly_with_idx(vector<int>& terminal, vector<int>& cand, \
// unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
// int kth, int& c2t){
//     int sink = terminal.back();
//     vector<vector<pair<int, int>>> disG;
//     int cnt = 0;
//     for(auto ele:umpii){
//         //neglect unchosen cand
//         if(0 != ele.second){
//             cnt++;
//         }
//     }
//     //terminal size+valid cand num.
//     disG.resize(cnt+1);

//     //cand assembly points and the sink.
//     //dis from cand to cand∪sink 
//     int cur = 0;
//     for(int l=0; l<kth; l++){
//         int loopCand = 0;
//         if(0 != umpii[cand[l]]){
//             //cand to cand
//             for(int m=0; m<kth; m++){
//                 if(0 != umpii[cand[m]]){
//                     disG[cur].push_back(make_pair(loopCand, distanceQuery(cand[loopCand], cand[cur])));
//                     loopCand++;
//                 }
//             }
//             //cand to sink
//             disG[cur].push_back(make_pair(loopCand, distanceQuery(cand[cur], sink)));
//             cur++;
//         }
//     }
//     //sink to cand
//     int candCnt = 0;
//     for(int i=0; i<kth; i++){
//         if(0!=umpii[cand[i]]){
//             disG[cur].push_back(make_pair(candCnt, distanceQuery(cand[candCnt], sink)));
//             candCnt++;
//         }
//     }
//     //sink to sink
//     disG[cur].push_back(make_pair(cur, distanceQuery(sink, sink)));
    
//     c2t = prim(disG, disG.size());
//     // cout<<"c2t: "<<c2t<<endl;
// }
// //algo2 
// void algo2_ComputePathForSource(vector<int>& terminal, vector<int>& cand, \
// unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
// int kth, int& s2c){
//     int sink = terminal.back();
//     //get dis from each term to the rest
//     for(int i=0; i<terminal.size(); i++){
//         dij(termDist[i], pre, npath, terminal[i]);
//     } 
//     //pass to algo3 and negelect the unchosen cand
//     for(int i=0; i<kth; i++){
//         umpii.emplace(make_pair(cand[i], 0));
//     }
//     //allocate to the nearest cand
//     for(int t=0; t<terminal.size()-1; t++){
//         int tmpSrc = terminal[t];
//         int src2sink = termDist[t][sink];
//         int sc2cd = INF;
//         int chosenCand = 0;
//         for(int i=0; i<kth; i++){
//             if(sc2cd > termDist[t][cand[i]]){
//                 sc2cd = termDist[t][cand[i]];
//                 chosenCand = i;
//             }
//         }
//         if(sc2cd>src2sink){
//             s2c+=src2sink;
//             // umpii[sink] += 1;
//         }else{
//             s2c+=sc2cd;
//             umpii[cand[chosenCand]] += 1;
//         }
//     }
//     cout<<"s2c:"<<s2c<<endl;
// }
// //algo3
// //cand without src will be neglected.
// void algo3_computePathAssembly(vector<int>& terminal, vector<int>& cand, \
// unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
// int kth, int& c2t){
//     vector<vector<pair<int, int>>> disG;
//     int cnt = 0;
//     for(auto ele:umpii){
//         //neglect unchosen cand
//         if(0 != ele.second){
//             cnt++;
//         }
//     }
//     //terminal size+valid cand num.
//     disG.resize(terminal.size()+cnt);

//     //terminal to terminal+cand
//     //t 0~terminal.size()-1
//     //cand newStart+cnt
//     int newStart = terminal.size();
//     for(int i=0; i<terminal.size(); i++){
//         //t2t
//         for(int j=0; j<terminal.size(); j++){
//             disG[i].push_back(make_pair(j,termDist[i][terminal[j]]));
//         }
//         int newCur = 0;
//         for(int l=0; l<kth; l++){
//             //neglect cand without src and can
//             // if(0 != umpii[cand[l]] && terminal[i] != cand[l]){
//             if(0 != umpii[cand[l]] ){
//                 disG[i].push_back(make_pair(newStart+newCur, termDist[i][cand[l]]));
//                 newCur++;
//             }
//         }
//     }

//     //cand to term and cand
//     vector<vector<int>> candDist(kth+1, vector<int>(GraphSize+1, INF));
//     vector<vector<int>> tmpPrev;
//     tmpPrev.resize(GraphSize+1);
//     vector<int> tmpnpath(GraphSize+1, 0);
//     for(int i=0; i<kth; i++){
//         if(0 != umpii[cand[i]]){
//             dij(candDist[i], tmpPrev, tmpnpath, cand[i]);
//         }
//     }
//     int newCur = 0;
//     for(int i=0; i<kth; i++){
//         if(0 != umpii[cand[i]]){
//             for(int t=0; t<terminal.size(); t++){
//                 //terminal to cand
//                 // if(cand[i] == terminal[t]) continue;
//                 if(0 == umpii[cand[i]]) continue;
//                 disG[newStart+newCur].push_back(make_pair(t, termDist[t][cand[i]]));
//             }
//             int tmpNewCur = 0;
//             for(int j=0; j<kth; j++){
//                 if(0 == umpii[cand[j]]) continue;
//                 disG[newStart+newCur].push_back(make_pair(newStart+tmpNewCur, candDist[i][cand[j]]));
//                 tmpNewCur++;
//             }
//             newCur++;
//         }
//     }
//     c2t = prim(disG, disG.size());
//     cout<<"c2t: "<<c2t<<endl;
// }

// //Algo1 is designed based on the fact we will exhuastively enumerate all possible combinations in the search tree.
// //Therefore, algo1 is not compatible with appx approach.
// //We propose a new greedy algo to greedily search the min cand in the search tree.
// void algo1_ComputeAllPath(vector<int>& terminal, vector<int>& cand, int budget, int idx_flag){
//     if(idx_flag){
//         int prevMinKth = INF;
//         //the cands are chosen already. 
//         //Therefore, no need to enum possible search procdure.
//         for (int kth = budget; kth <= budget; kth++)
//         {
//             //pair<candId, assigned srcs>
//             unordered_map<int, int> umpii;
//             vector<vector<int>> termDist;
//             // path loaded but not retrieved.
//             vector<vector<int>> pre;
//             pre.resize(GraphSize + 1);
//             vector<int> npath(GraphSize + 1, 0);
//             auto kth_start = timerSt(); 
//             int s2c = 0;
//             int c2t = 0;
//             algo2_ComputePathForSource_with_idx(terminal, cand, umpii, termDist, pre, npath, kth, s2c); 
//             if(s2c >= prevMinKth) continue;
//             algo3_computePathAssembly_with_idx(terminal, cand, umpii, termDist, pre, npath, kth, c2t);
//             if(prevMinKth > s2c+c2t){
//                 prevMinKth = s2c+c2t;
//             }
//             cout<<"kth: "<<kth<<" dist: "<<prevMinKth<<endl;
//             cout<<"cands and if chosen: "<<endl;
//             for(auto ele:umpii){
//                 cout<<ele.first<<" "<<ele.second<<"| ";  
//             }
//             cout<<endl;
//             auto kth_end = timerNd(); 
//             double dur_time = dur(kth_end, kth_start);
//             res50Time[kth].push_back(dur_time);
//             res50Dist[kth].push_back(prevMinKth);
//             // cout<<"k="<<kth<<" dist "<<(s2c+c2t)<<endl;
//         }
//     }else{
//         //cands are chosen already.
//         for (int kth = budget; kth <= budget; kth++)
//         {
//             // cout << "*_*_*_*" << kth << "*_*_*_*" << endl;
//             unordered_map<int, int> umpii;
//             vector<vector<int>> termDist(terminal.size() + 1, vector<int>(GraphSize + 1, INF));
//             // path loaded but not retrieved.
//             vector<vector<int>> pre;
//             pre.resize(GraphSize + 1);
//             vector<int> npath(GraphSize + 1, 0);
//             clock_t kth_start = clock();
//             int s2c = 0;
//             int c2t = 0;
//             algo2_ComputePathForSource(terminal, cand, umpii, termDist, pre, npath, kth, s2c);
//             algo3_computePathAssembly(terminal, cand, umpii, termDist, pre, npath, kth, c2t);
//             clock_t kth_end = clock();
//             clock_t tot_kth = kth_end - kth_start;
//             res50Time[kth].push_back((double)(1.0 * tot_kth / CLOCKS_PER_SEC));
//             res50Dist[kth].push_back(s2c + c2t);
//         }
//         // clock_t algo2_nd = clock();
//         // cout << "algo 2: " << (double)(algo2_nd - algo2_st) / CLOCKS_PER_SEC << endl;
//     }
// }
// int algo_simpleGreedy(vector<int>& terminal, int budget, vector<int>& aps, queue<int>& qi, int& upDist){
//     vector<int> src(terminal.begin(), terminal.end()-1);
//     int sink = terminal.back();
//     vector<vector<int>> vvi_times(terminal.size(), vector<int>(GraphSize+1, 0));
//     vector<int> totDist(GraphSize+1, 0);
//     vector<int> longest_edge_of_src(terminal.size()+1, 0);
//     int cliqueSize = terminal.size()+budget;
//     vector<vector<pair<int, int>>> newG(cliqueSize, vector<pair<int, int>>(cliqueSize, make_pair(0, 0)));

//     if(budget == 1){
//         // queue<int> qi;
//         int minDist = 0;
//         // tot dist should be no larger than direct to sink
//         for (auto ele : src)
//         {
//             minDist += distanceQuery(ele, sink);
//         }
//         cout << "upperbound: " << minDist << endl;
//         int first_ap = 0;
//         int first_ap_cnt = 0;
//         for (int i = 1; i < GraphSize + 1; i++)
//         {
//             if (i == sink)
//                 continue;
//             int tmp_ap_cnt = 0;
//             int curDist = distanceQuery(i, sink);
//             for (int j = 0; j < src.size(); j++)
//             {
//                 int src2i = distanceQuery(i, src[j]);
//                 int src2t = distanceQuery(src[j], sink); 
//                 if (src2i > src2t)
//                 {
//                     curDist += src2t;
//                 }
//                 else
//                 {
//                     curDist += src2i;
//                     tmp_ap_cnt++;
//                 }
//             }
//             if (minDist > curDist)
//             {
//                 first_ap = i;
//                 minDist = curDist;
//                 first_ap_cnt = tmp_ap_cnt;
//                 // cout << "curDist: " << curDist << endl;
//                 qi.push(i);
//             }
//         }
//         cout << "minDist 1st k: " << minDist << endl;
//         for (int i = 0; i < src.size(); i++)
//         {
//             longest_edge_of_src[i] = min(distanceQuery(sink, src[i]),\
//             distanceQuery(first_ap, src[i]));
//         }
//         int fir_par = 0;
//         for(auto ele:longest_edge_of_src){
//             fir_par+=ele;
//         }
//         cout<<"fir par weight: "<<fir_par<<endl;
//         aps.push_back(first_ap);
//         cout << "1st ap cnt: " << first_ap_cnt << endl;

//         if (budget <= 1 || first_ap_cnt == 0)
//         {
//             return minDist;
//         }
//     }
//     //k>1
//     int qisize = qi.size();
//     int cnt = 0;
//     int kth = 1;
//     int kthMinDist = upDist;
//     int kthAP = 0;
//     while(!qi.empty()){
//         int u = qi.front();
//         qi.pop();
//         cnt++;
//         if(cnt == qisize){
//             //no valid next ap, break
//             if(kthMinDist >= upDist){
//                 break;
//             }
//             // upd upperbound...
//             upDist = kthMinDist;
//             aps.push_back(kthAP);
//             // kth ...
//             kth++;
//             if(kth >= budget){
//                 break;
//             }
//         }
//         int excess_cnt = 0;
//         for(int j=0; j<src.size(); j++){
//             int s2u = distanceQuery(src[j], u);
//             if(s2u > longest_edge_of_src[j]){
//                 excess_cnt++;
//             }
//         }
//         //simple pruning, remove unpromising cands.
//         if(excess_cnt < 4){
//             //possible cand for next ap
//             qi.push(u);
//             int s2c=0;
//             //assgin to nearest 
//             for(int i=0; i<src.size(); i++){
//                 int src2cand = distanceQuery(u, src[i]);
//                 //assigned to either u
//                 if(src2cand < longest_edge_of_src[i]){
//                     s2c += src2cand;
//                 //or stay unchanged.
//                 }else{
//                     s2c += longest_edge_of_src[i];
//                 }
//             }
//             cout<<"s2c: "<<s2c<<endl;
//             // make clique: src, aps, u
//             aps.push_back(u);
//             fillG(aps, sink, newG);
//             //prim
//             int c2t = prim(newG, aps.size()+1);
//             aps.pop_back();
//             //comp if valid
//             if(c2t+s2c < kthMinDist){
//                 kthMinDist = c2t+s2c;
//                 kthAP = u;
//             }
//         }
//     }

//     return upDist;
// }
// // loop over all vertices and find next possible ap. 
// // O(kn), it is ugly slow.
// int algo_greedy(vector<int>& terminal, vector<int>& aps, int budget){
//     // init
//     vector<int> src(terminal.begin(), terminal.end()-1);
//     int cliqueSize = terminal.size()+budget;
//     vector<vector<pair<int, int>>> newG(cliqueSize, vector<pair<int, int>>(cliqueSize, make_pair(0, 0)));
//     int sink = terminal.back();
//     int minDist = 0;
//     for(auto ele:src){
//         minDist += distanceQuery(ele, sink);
//     }
//     int ap=0;
//     int loopMinDist = minDist;
//     while(1){
//         int curDist = 0;
//         for(int i=1; i<GraphSize+1; i++){
//             if(i%10000 == 0){
//                 cout<<i<<endl;
//             }
//             if (i == sink)
//                 continue;
//             int curDist = 0; 
//             aps.push_back(i);
//             for (int j = 0; j < src.size(); j++)
//             {
//                 int src2t = distanceQuery(src[j], sink);
//                 int src2i = INF; 
//                 int whichaps = 0;
//                 for(int kth=0; kth<aps.size(); kth++){
//                     int src2loopc =  distanceQuery(aps[kth], src[j]);
//                     if(src2i > src2loopc){
//                         src2i = src2loopc;
//                     }
//                 }
//                 if (src2i > src2t)
//                 {
//                     curDist += src2t;
//                 }
//                 else
//                 {
//                     curDist += src2i;
//                 }
//             }
//             //if 1st part is greater than upbound, no further computation for 2nd part.
//             if(curDist > loopMinDist) continue;
//             cout<<"src+aps: "<<src.size()+aps.size()<<endl;
//             fillG(aps, sink, newG);
//             int c2t = prim(newG, src.size()+aps.size());
//             curDist += c2t;
//             if(loopMinDist > curDist){
//                 ap = i;
//                 loopMinDist = curDist;
//             }
//             aps.pop_back();
//         }//end for
//         if(minDist > loopMinDist){
//             aps.push_back(ap);
//             minDist = loopMinDist;
//         }else{
//             break; 
//         }
//         cout<<"minDist: "<<minDist<<" kth: "<<aps.size()<<endl;
//     }//end while
//     return minDist;
// }
// //worst case O(n^2)
// //It runs slowly for for larger k
// void algo_apx_greedy(vector<int>& terminal, vector<int>& cand, int budget){
//     res50Dist.resize(budget+1);
//     res50Time.resize(budget+1);
//     int sink = terminal.back();
//     vector<int> src(terminal.begin(), terminal.end()-1);
//     long long ub = 0;
//     for(auto ele:src){
//         ub+=distanceQuery(ele, sink);
//     }
//     cout<<"ub: "<<ub<<endl;
//     int pos = 0;
//     int newStartPos = 0;
//     long long minDist = LLINF;
//     long long curMinDist = LLINF;
//     vector<int> aps;
//     //in the search tree, we greedily select the most promising cand each time.
//     //Such that we do not need to enum all the possible combinations as in BnB.
//     //i starts from 0 because cand starts from zero
//     for(int i=0; i<budget; i++){
//         auto apx_grdy_start = timerSt();
//         for(int j=newStartPos; j<budget; j++){
//             aps.push_back(cand[j]);
//             vector<vector<pair<int, int>>> newG;
//             // auto s2c_st = timerSt();
//             unordered_map<int, vector<int>> umpii;
//             long long s2c = assign_src2cand(terminal, aps, newG, umpii);
//             // auto s2c_nd = timerNd();
//             // auto c2t_st = timerSt();
//             long long c2t = prim(newG, newG.size());
//             // auto c2t_nd = timerNd();
//             // double s2c_dur = dur(s2c_nd, s2c_st);
//             // double c2t_dur = dur(c2t_nd, c2t_st);
//             // cout<<"s2c_dur: "<<s2c_dur<<" c2t_dur: "<<c2t_dur<<endl;
//             if(curMinDist > s2c+c2t){
//                 curMinDist = s2c+c2t;
//                 newStartPos = j;
//             }
//             aps.pop_back();
//         }
//         // cout<<"curMinDist: "<<curMinDist<<" minDist: "<<minDist<<endl;
//         //collect res.
//         if(minDist >= curMinDist){
//             minDist = curMinDist;
//             aps.push_back(cand[newStartPos]);
//         }
//         if(minDist < curMinDist && i>0||newStartPos+1>=budget){
//             for(int st=i; st<budget; st++){
//                 res50Dist[st+1].push_back(minDist);
//                 res50Time[st+1].push_back(1.0*MYMAX_INT_HALF);
//             }
//             break;
//         }
//         //next unchosen cand
//         newStartPos++; 
//         curMinDist = INF;
//         auto apx_grdy_end = timerNd();
//         double dur_time = dur(apx_grdy_end, apx_grdy_start);
//         res50Dist[i+1].push_back(minDist);
//         res50Time[i+1].push_back(dur_time);
//     }
//     cout<<"minDist apx: "<<minDist<<endl
//     <<"ap size: "<<aps.size()<<endl
//     <<"aps: "<<endl;
//     for(auto ele:aps){
//         cout<<ele<<" ";
//     }
//     cout<<endl;
// }

