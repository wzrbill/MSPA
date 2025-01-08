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
double tot_dij_time = 0.0;
double tot_idx_time = 0.0;
double tot_get_paths_time = 0.0;
double tot_other_1 = 0.0;
double tot_other_2 = 0.0;
double tot_path_time = 0.0;
// int minWeight = 0;
vector<vector<long long>> res50Dist(6);
vector<vector<double>> res50Time(6);

// PrunedLandmarkLabeling<> pll;
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

//RECOVER_MSPATREE
struct Pent{
    int nodeid;
    int s1;
    int r1;
    int s2;
    int r2;
    int isAP;
    Pent():nodeid(0),s1(0),r1(0),s2(0),r2(0),isAP(0){}  
    Pent(int n, int sub, int rsub, int cmpSub, int rcmpSub, int ap):nodeid(n), s1(sub), r1(rsub), s2(cmpSub), r2(rcmpSub), isAP(ap){}
};
//RECOVER_MSAPATREE
vector<vector<vector<vector<Pent>> > > pre;
vector<vector<vector<vector<int> > > > foot;
vector<Pent> vp; 
vector<pair<int, int>> treeEdge;;

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
void dfsTrees(vector<int>& Srcs, int nodeid, int state, int r, int ap, int pos){
    if(nodeid == 0) return;
    if(r<ap) return;
    cout<<"["<<state<<" ,"<<nodeid<<", "<<r<<", "<<ap<<"] = "
        <<"["<<pre[state][nodeid][r][ap].s1
        <<", "<< pre[state][nodeid][r][ap].nodeid
        <<", "<< pre[state][nodeid][r][ap].r1
        <<", "<<pre[state][nodeid][r][ap].isAP<<"] + "
        <<" ["<<pre[state][nodeid][r][ap].s2
        <<", "<< pre[state][nodeid][r][ap].nodeid
        <<", "<< pre[state][nodeid][r][ap].r2
        <<", "<<pre[state][nodeid][r][ap].isAP<<"] "<<endl;

    if(foot[state][nodeid][r][ap] == 1) return;
    /* if(preState[state][nodeid].nodeid == -1) return; */
    foot[state][nodeid][r][ap] = 1;       
    int preNode = pre[state][nodeid][r][ap].nodeid;
    int s1 = pre[state][nodeid][r][ap].s1;
    int r1 = pre[state][nodeid][r][ap].r1;
    int s2 = pre[state][nodeid][r][ap].s2;
    int r2 = pre[state][nodeid][r][ap].r2;
    int isap = pre[state][nodeid][r][ap].isAP;
    /* cout<<"pre "<<preNode<<" s1 "<<s1<<" s2 "<<s2<<endl; */
    if(nodeid == preNode){
        dfsTrees(Srcs, preNode, s1, r1, isap, 2*pos+1);
        vp[2*pos+1] = pre[s1][preNode][r1][isap];
        dfsTrees(Srcs, preNode, s2, r2, isap, 2*pos+2);
        vp[2*pos+2] = pre[s2][preNode][r2][isap];
    }else{
        dfsTrees(Srcs, preNode, s1, r1, isap, 2*pos+1);
        if(preNode == 0){
            treeEdge.push_back(make_pair(nodeid, Srcs[(int)log2(state)]));
        }else{
            treeEdge.push_back(make_pair(nodeid, preNode));
        }
        vp[2*pos+1] = pre[s1][preNode][r1][isap];
    }
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

void printRes2File(){
    for(int i=1; i<(int)res50Dist.size(); i++){
        if(res50Dist[i].size() == 0) {
            // cout<<i<<"th res is void"<<endl;
            continue;
        }
        // cout<<i<<"th size "<<res50Dist[i].size()<<endl; 
        for(int j=0; j<(int)res50Dist[i].size(); j++){
            fprintf(outOpenFile, "%d", res50Dist[i][j]); 
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
int timeOutCheck(std::chrono::time_point<std::chrono::high_resolution_clock> algo_st){
    int flag = 0;
    auto checkTimer = timerNd();
    double timeout = dur(checkTimer, algo_st);
    if (timeout >= 7200.0){
        flag = 1;
    }
    return flag;
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
long long terminals_cand_tot_dis(vector<int>& terminal, int cand){
    long long t_tot_dis = 0;
    for(auto ele:terminal){
        t_tot_dis += distanceQuery(ele, cand);
    }
    return t_tot_dis;
}
int dij(int src, int sink){
    vector<int> termDist(GraphSize+1, INF);
    termDist[src]=0;
    vector<int> vis(GraphSize+1, 0);
    priority_queue<pair<int, int> > pq_local;
    pq_local.emplace(make_pair(0, src));
    while(!pq_local.empty()){
        pair<int, int> pii = pq_local.top();
        int u = pii.second;
        pq_local.pop();
        if(vis[u]) continue;
        vis[u] = 1;
        for(int kth=0; kth<Degree[u]; kth++){
            int nei = Neighbor[u][kth];
            int newDist = termDist[u] + Weight[u][kth];
            if(termDist[nei] > newDist){
                termDist[nei] = newDist;
                pq_local.emplace(make_pair(newDist, nei));
            }
        }
    }
    return termDist[sink];
}
void dij( vector<int>& termDist, vector<vector<int>>& pre, vector<int>& npath, int src){
// void dij( vector<int>& termDist,  vector<int>& npath, int src){
    // init distances
    npath[src] = 1;
    termDist[src] = 0;
    auto initPath = [&](int u, int v) {
        if(0 != pre.size() ){
            pre[v] = {u};
        }
            npath[v] = npath[u];
    };
    vector<int> vis(GraphSize+1, 0);
    priority_queue<pair<int, int> > pq_local;
    pq_local.emplace(make_pair(0, src));
    while(!pq_local.empty()){
        pair<int, int> pii = pq_local.top();
        int u = pii.second;
        pq_local.pop();
        if(vis[u]) continue;
        vis[u] = 1;
        for(int kth=0; kth<Degree[u]; kth++){
            int nei = Neighbor[u][kth];
            int newDist = termDist[u] + Weight[u][kth];
            if(termDist[nei] > newDist){
                initPath(u, nei);
                termDist[nei] = newDist;
                pq_local.emplace(make_pair(newDist, nei));
            }else if(termDist[nei] == newDist){
                if(0 != pre.size()){
                    pre[nei].push_back(u);
                }
                npath[nei] += npath[u];
            }
        }
    }

}