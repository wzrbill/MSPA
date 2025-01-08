#include "MSPAUtilities.h"
#include "pruned_landmark_labeling.h"
#include "H2H.h"
#include "appx.h"
using namespace std;
// x Deprecated
//algo2 with idx, no path retrieval.
void algo2_ComputePathForSource_with_idx(vector<int>& terminal, vector<int>& cand, \
unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
int kth, int& s2c){
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    //will be delivered to algo3, and will negelect the unchosen cand
    for(int i=0; i<kth; i++){
        umpii.emplace(make_pair(cand[i], 0));
    }
    //allocate to the nearest cand
    for(int t=0; t<src.size(); t++){
        int tmpSrc = src[t];
        int src2sink = distanceQuery(t, sink);
        int sc2cd = INF;
        int chosenCand = 0;
        //choose cand limited by budget
        for(int i=0; i<kth; i++){
            int tmpSrc2cand = distanceQuery(tmpSrc, cand[i]);
            if(sc2cd > tmpSrc2cand){
                sc2cd = tmpSrc2cand; 
                chosenCand = i;
            }
        }
        if(sc2cd>src2sink){
            s2c+=src2sink;
        // if not assigned to sink, then mark cand chosen by one src.
        }else{
            s2c+=sc2cd;
            umpii[cand[chosenCand]] += 1;
        }
    }
    // cout<<"s2c:"<<s2c<<endl;
}
// x Deprecated
//algo3, no path retrieval.
void algo3_computePathAssembly_with_idx(vector<int>& terminal, vector<int>& cand, \
unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
int kth, int& c2t){
    int sink = terminal.back();
    vector<vector<pair<int, int>>> disG;
    int cnt = 0;
    for(auto ele:umpii){
        //neglect unchosen cand
        if(0 != ele.second){
            cnt++;
        }
    }
    //terminal size+valid cand num.
    disG.resize(cnt+1);

    //cand assembly points and the sink.
    //dis from cand to cand∪sink 
    int cur = 0;
    for(int l=0; l<kth; l++){
        int loopCand = 0;
        if(0 != umpii[cand[l]]){
            //cand to cand
            for(int m=0; m<kth; m++){
                if(0 != umpii[cand[m]]){
                    disG[cur].push_back(make_pair(loopCand, distanceQuery(cand[loopCand], cand[cur])));
                    loopCand++;
                }
            }
            //cand to sink
            disG[cur].push_back(make_pair(loopCand, distanceQuery(cand[cur], sink)));
            cur++;
        }
    }
    //sink to cand
    int candCnt = 0;
    for(int i=0; i<kth; i++){
        if(0!=umpii[cand[i]]){
            disG[cur].push_back(make_pair(candCnt, distanceQuery(cand[candCnt], sink)));
            candCnt++;
        }
    }
    //sink to sink
    disG[cur].push_back(make_pair(cur, distanceQuery(sink, sink)));
    
    c2t = prim(disG, disG.size());
    // cout<<"c2t: "<<c2t<<endl;
}
// x Deprecated
//algo2 
//todo allow retrieve paths: getPath(). Here we only have path stored but not retrieved.
void algo2_ComputePathForSource(vector<int>& terminal, vector<int>& cand, \
unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
int kth, int& s2c){
    int sink = terminal.back();
    //get dis from each term to the rest
    for(int i=0; i<terminal.size(); i++){
        dij(termDist[i], pre, npath, terminal[i]);
    } 
    //pass to algo3 and negelect the unchosen cand
    for(int i=0; i<kth; i++){
        umpii.emplace(make_pair(cand[i], 0));
    }
    //allocate to the nearest cand
    for(int t=0; t<terminal.size()-1; t++){
        int tmpSrc = terminal[t];
        int src2sink = termDist[t][sink];
        int sc2cd = INF;
        int chosenCand = 0;
        for(int i=0; i<kth; i++){
            if(sc2cd > termDist[t][cand[i]]){
                sc2cd = termDist[t][cand[i]];
                chosenCand = i;
            }
        }
        if(sc2cd>src2sink){
            s2c+=src2sink;
            // umpii[sink] += 1;
        }else{
            s2c+=sc2cd;
            umpii[cand[chosenCand]] += 1;
        }
    }
    cout<<"s2c:"<<s2c<<endl;
}
//x Deprecated
//algo3
//cand without src will be neglected.
//TODO allow retrieve paths. Use pre to recover path.
void algo3_computePathAssembly(vector<int>& terminal, vector<int>& cand, \
unordered_map<int, int>& umpii, vector<vector<int>>& termDist, vector<vector<int>>& pre, vector<int>& npath,\
int kth, int& c2t){
    vector<vector<pair<int, int>>> disG;
    int cnt = 0;
    for(auto ele:umpii){
        //neglect unchosen cand
        if(0 != ele.second){
            cnt++;
        }
    }
    //terminal size+valid cand num.
    disG.resize(terminal.size()+cnt);

    //terminal to terminal+cand
    //t 0~terminal.size()-1
    //cand newStart+cnt
    int newStart = terminal.size();
    for(int i=0; i<terminal.size(); i++){
        //t2t
        for(int j=0; j<terminal.size(); j++){
            disG[i].push_back(make_pair(j,termDist[i][terminal[j]]));
        }
        int newCur = 0;
        for(int l=0; l<kth; l++){
            //neglect cand without src and can
            // if(0 != umpii[cand[l]] && terminal[i] != cand[l]){
            if(0 != umpii[cand[l]] ){
                disG[i].push_back(make_pair(newStart+newCur, termDist[i][cand[l]]));
                newCur++;
            }
        }
    }

    //cand to term and cand
    vector<vector<int>> candDist(kth+1, vector<int>(GraphSize+1, INF));
    vector<vector<int>> tmpPrev;
    tmpPrev.resize(GraphSize+1);
    vector<int> tmpnpath(GraphSize+1, 0);
    for(int i=0; i<kth; i++){
        if(0 != umpii[cand[i]]){
            dij(candDist[i], tmpPrev, tmpnpath, cand[i]);
        }
    }
    int newCur = 0;
    for(int i=0; i<kth; i++){
        if(0 != umpii[cand[i]]){
            for(int t=0; t<terminal.size(); t++){
                //terminal to cand
                // if(cand[i] == terminal[t]) continue;
                if(0 == umpii[cand[i]]) continue;
                disG[newStart+newCur].push_back(make_pair(t, termDist[t][cand[i]]));
            }
            int tmpNewCur = 0;
            for(int j=0; j<kth; j++){
                if(0 == umpii[cand[j]]) continue;
                disG[newStart+newCur].push_back(make_pair(newStart+tmpNewCur, candDist[i][cand[j]]));
                tmpNewCur++;
            }
            newCur++;
        }
    }
    c2t = prim(disG, disG.size());
    cout<<"c2t: "<<c2t<<endl;
}

//x Deprecated
//Algo1 is designed based on the fact we will exhuastively enumerate all possible combinations in the search tree.
//Therefore, algo1 is not compatible with appx approach.
//We propose a new greedy algo to greedily search the min cand in the search tree.
void algo1_ComputeAllPath(vector<int>& terminal, vector<int>& cand, int budget, int idx_flag){
    if(idx_flag){
        int prevMinKth = INF;
        //the cands are chosen already. 
        //Therefore, no need to enum possible search procdure.
        for (int kth = budget; kth <= budget; kth++)
        {
            //pair<candId, assigned srcs>
            unordered_map<int, int> umpii;
            vector<vector<int>> termDist;
            // path loaded but not retrieved.
            vector<vector<int>> pre;
            pre.resize(GraphSize + 1);
            vector<int> npath(GraphSize + 1, 0);
            auto kth_start = timerSt(); 
            int s2c = 0;
            int c2t = 0;
            algo2_ComputePathForSource_with_idx(terminal, cand, umpii, termDist, pre, npath, kth, s2c); 
            if(s2c >= prevMinKth) continue;
            algo3_computePathAssembly_with_idx(terminal, cand, umpii, termDist, pre, npath, kth, c2t);
            if(prevMinKth > s2c+c2t){
                prevMinKth = s2c+c2t;
            }
            cout<<"kth: "<<kth<<" dist: "<<prevMinKth<<endl;
            cout<<"cands and if chosen: "<<endl;
            for(auto ele:umpii){
                cout<<ele.first<<" "<<ele.second<<"| ";  
            }
            cout<<endl;
            auto kth_end = timerNd(); 
            double dur_time = dur(kth_end, kth_start);
            res50Time[kth].push_back(dur_time);
            res50Dist[kth].push_back(prevMinKth);
            // cout<<"k="<<kth<<" dist "<<(s2c+c2t)<<endl;
        }
    }else{
        //cands are chosen already.
        //idx_flag is 0 and no index
        for (int kth = budget; kth <= budget; kth++)
        {
            // cout << "*_*_*_*" << kth << "*_*_*_*" << endl;
            unordered_map<int, int> umpii;
            vector<vector<int>> termDist(terminal.size() + 1, vector<int>(GraphSize + 1, INF));
            // path loaded but not retrieved.
            vector<vector<int>> pre;
            pre.resize(GraphSize + 1);
            vector<int> npath(GraphSize + 1, 0);
            clock_t kth_start = clock();
            int s2c = 0;
            int c2t = 0;
            algo2_ComputePathForSource(terminal, cand, umpii, termDist, pre, npath, kth, s2c);
            algo3_computePathAssembly(terminal, cand, umpii, termDist, pre, npath, kth, c2t);
            clock_t kth_end = clock();
            clock_t tot_kth = kth_end - kth_start;
            res50Time[kth].push_back((double)(1.0 * tot_kth / CLOCKS_PER_SEC));
            res50Dist[kth].push_back(s2c + c2t);
        }
        // clock_t algo2_nd = clock();
        // cout << "algo 2: " << (double)(algo2_nd - algo2_st) / CLOCKS_PER_SEC << endl;
    }
}


//x Deprecated
//we traverse the search tree with reordered ver dis to all srcs in G'
//the upperbound is only valid for its subtrees () not the cousin trees
void search(vector<int>& apset, int k, vector<pair<int, long long>>& Gprim, vector<int>& terminal, int pos, int upBd, vector<int>& rstar){
    int given_ap_tot_dis = algo1_comp_all_path(apset, terminal, Gprim, 1);
    //when apset size is 0, we use upperbound as cur tot dis
    if(-1 == given_ap_tot_dis){
        given_ap_tot_dis = upBd;
    }
    //invalid ap set return
    if(given_ap_tot_dis > upBd){
        return;
    }
    if(apset.size()<k){
        for(int i=pos; i<GraphSize; i++){
            apset.push_back(Gprim[i].first);
            //if it is bounded, keep searching
            if(lowerBound(apset, terminal) < given_ap_tot_dis){
                //cur ver is regarded as a root to the current search tree, so cur tot dis is upBd
                search(apset, k, Gprim, terminal, i+1, given_ap_tot_dis, rstar);
            }
            apset.pop_back();
        }
    }else if(apset.size() == k){
        if(given_ap_tot_dis < upBd){
            cout<<upBd<<endl;
            if(0 == rstar.size()){
                rstar = apset;
            }else{
                if(given_ap_tot_dis < algo1_comp_all_path(rstar, terminal, Gprim, 1)){
                    rstar = apset;
                }
            }
            for(auto ele:rstar){
                cout<<ele<<" ";
            }
            cout<<endl;
        }
    }
}
//x Deprecated
void algo4_Branch_and_Bound(vector<int>& terminals, int k){
    //reorder by dis
    vector<pair<int, long long>> reorder_by_dis;
    for(int i=1; i<=GraphSize; i++){
        reorder_by_dis.push_back({i, terminals_cand_tot_dis(terminals, i)});
    }
    sort(reorder_by_dis.begin(), reorder_by_dis.end(), [&](pair<int, int> pa, pair<int, int> pb){
        return pa.second<pb.second;
    });

    //tmp set of assembly points during comp
    vector<int> apset;
    int upBd = 0;
    int sink = terminals.back();
    vector<int> srcs(terminals.begin(), terminals.end()-1);
    //any dis longer than direct tot dis to des will be pruned
    for(auto ele:srcs){
        upBd += distanceQuery(ele, sink);
    }
    //optimal solution of assembly points
    vector<int> rstar;

    //all init done, start search on the search tree
    search(apset, k, reorder_by_dis, terminals, 0, upBd, rstar);
}

int algo1_comp_all_path(vector<int>& apset, vector<int>& terminal, vector<pair<int, long long>>& Gprim, int code){
    //no cand in apset, we skip this computation
    if(0 == apset.size()){
        return -1;
    }
    int tot_dis=0;
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    int s2c = 0;
    int c2t = 0;
    //algo2
    //we assign each src to its nearest assembly point or the des
    for(int i=0; i<srcs.size(); i++){
        int min_Src2cand = code==1?distanceQuery(srcs[i], sink):dij(srcs[i], sink); 
        for(int j=0; j<apset.size(); j++){
            int src2cand = code==1?distanceQuery(apset[j], srcs[i]):dij(srcs[i], sink);
            if(min_Src2cand > src2cand){
                min_Src2cand = src2cand;
            }
        }
        s2c+=min_Src2cand;
    }
    apset.push_back(sink);
    vector<vector<pair<int, int>>> newG;
    //all pairs distance matrix
    fillG(newG, apset);
    //algo3
    c2t = prim(newG, apset.size());
    apset.pop_back();
    tot_dis += s2c;
    tot_dis += c2t;
    return tot_dis;
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
void search_tree(vector<int>& apset, int k, vector<pair<int, long long>>& Gprim, vector<int>& terminal, int pos, int upBd, vector<int>& rstar, int code){
    int given_ap_tot_dis = algo1_comp_all_path(apset, terminal, Gprim, code);
    //when apset size is 0, we use upperbound as cur tot dis
    if(-1 == given_ap_tot_dis){
        given_ap_tot_dis = upBd;
    }
    //invalid ap set return
    if(given_ap_tot_dis > upBd){
        return;
    }
    if(apset.size()<k){
        for(int i=pos; i<GraphSize; i++){
            apset.push_back(Gprim[i].first);
            //if it is bounded, keep searching
            if(lowerBound(apset, terminal) < given_ap_tot_dis){
                //cur ver is regarded as a root to the current search tree, so cur tot dis is upBd
                search(apset, k, Gprim, terminal, i+1, given_ap_tot_dis, rstar);
            }
            apset.pop_back();
        }
    }else if(apset.size() == k){
        if(given_ap_tot_dis < upBd){
            cout<<upBd<<endl;
            if(0 == rstar.size()){
                rstar = apset;
            }else{
                if(given_ap_tot_dis < algo1_comp_all_path(rstar, terminal, Gprim, code)){
                    rstar = apset;
                }
            }
            for(auto ele:rstar){
                cout<<ele<<" ";
            }
            cout<<endl;
        }
    }
}


void algo4_BnB(vector<int>& terminals, int k, int code){
    //reorder by dis
    vector<pair<int, long long>> reorder_by_dis;
    for(int i=1; i<=GraphSize; i++){
        reorder_by_dis.push_back({i, terminals_cand_tot_dis(terminals, i)});
    }
    sort(reorder_by_dis.begin(), reorder_by_dis.end(), [&](pair<int, int> pa, pair<int, int> pb){
        return pa.second<pb.second;
    });

    //tmp set of assembly points during comp
    vector<int> apset;
    int upBd = 0;
    int sink = terminals.back();
    vector<int> srcs(terminals.begin(), terminals.end()-1);
    //any dis longer than direct tot dis to des will be pruned
    for(auto ele:srcs){
        upBd += distanceQuery(ele, sink);
    }
    //optimal solution of assembly points
    vector<int> rstar;

    //all init done, start search on the search tree
    search_tree(apset, k, reorder_by_dis, terminals, 0, upBd, rstar, code);
}

void algo5_BasicDP(vector<int>& terminals, int k){
    int minDist = INF;
    //num of terminals
    int t = terminals.size();
    //num of src
    int sNum = t-1;
    //full mask
    int full_mask = (1<<sNum)-1;
    //nums of nodes;
    int v = GraphSize;
    vector<int> src(terminals.begin(), terminals.end()-1);
    int sink = terminals.back();
    cout<<"srcs: ";
    for(auto ele:src){
        cout<<ele<<" ";
    }
    cout<<" sink "<<sink<<endl;
    // cout<<"full_mask "<<full_mask<<" bin "<<hex<<full_mask<<dec<<endl;
    //[mask][nodid][rpCnt]
    //mask 1~2^(t-1)-1 ; nodeid 1~n; rpCnt 0~r;
    // vector<vector<vector<int> > > dp(1<<sNum, 
    //            vector<vector<int> > (v+1, 
    //                      vector<int>(bitCnt((1<<sNum)-1)+1, INF)));

    int bmsk = 1<<sNum;
    int vnum = v+1;
    int bgt = bitCnt((1<<sNum)-1)+1;
    // cpp style
    // int*** dp = new int**[bmsk];
    // for (int i = 0; i < bmsk; i++) {
    //     dp[i] = new int*[vnum];
    //     for (int j = 0; j < vnum; j++) {
    //         dp[i][j] = new int[bgt];
    //     }
    // }
    
    // // Don't forget to free the memory when done
    // for (int i = 0; i < bmsk; i++) {
    //     for (int j = 0; j < vnum; j++) {
    //         delete[] dp[i][j];
    //     }
    //     delete[] dp[i];
    // }
    // delete[] dp;
    int ***dp = (int ***)malloc(bmsk * sizeof(int **));
    for(int i=0; i<bmsk; i++){
        dp[i] = (int **)malloc(vnum * sizeof(int *));
        for(int j=0; j<vnum; j++){
            dp[i][j] = (int *)malloc(bgt * sizeof(int ));
        }
    }
    for(int i=0; i<bmsk; i++){
        for(int j=0;j<vnum; j++){
            for(int kth=0; kth<bgt; kth++){
                dp[i][j][kth] = INF;
            }
        }
    }
#ifdef RECOVER_MSPATREE
    prevTriState*** triTree = new prevTriState**[bmsk];
    for(int i=0; i<bmsk; i++){
        triTree[i] = new prevTriState*[vnum];
        for(int j=0; j<vnum; j++){
            triTree[i][j] = new prevTriState[bgt];
        }
    }
#endif
    vector<clock_t> timeCostEachAssemblyBudget(k+2, 0);

    //@caution nodeid starts from 1
    clock_t st = clock();
    for(int mask = 1; mask<1<<(sNum); mask++){
        dp[mask][sink][0] = bitDist(src, mask, sink);
        for(int i = 1; i<v+1; i++){
            if(i == sink) continue;
            dp[mask][i][1] = bitDist(src, mask, i);
#ifdef RECOVER_MSPATREE
            triTree[mask][i][1] = prevTriState(0, 0, 0, 0, 0);
#endif
        }
    }
 
    
    clock_t et = clock();
    timeCostEachAssemblyBudget[1] = et-st+1;
    //tree merge & tree grow  
    cout<<"init fin"<<endl;
    for(int r=2; r<=k; r++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
            for(int i=1;i<v+1; i++){
                if(i == sink) continue;
                // tree merge
                for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
                // for(int submask=1; submask<mask; submask++){
                    for(int subr=1; subr<=r; subr++){
                        if(0 == (r+1-subr)) continue;
                        if(dp[mask][i][r] > dp[submask][i][subr]+dp[submask^mask][i][(r+1)-subr]){
                            dp[mask][i][r] = dp[submask][i][subr]+dp[submask^mask][i][(r+1)-subr];
#ifdef RECOVER_MSPATREE
                            triTree[mask][i][r] = prevTriState(i, subr, submask, (r+1)-subr, submask^mask);
#endif
                        }
                    }
                }
                // tree grow
                for(int u=1; u<v+1; u++){
                    if(dp[mask][i][r] > dp[mask][u][r-1]+distanceQuery(u, i)){
                        dp[mask][i][r] = dp[mask][u][r-1]+distanceQuery(u, i);
#ifdef RECOVER_MSPATREE
                            triTree[mask][i][r] = prevTriState(u, r-1, mask, 0, 0);
#endif
                    }
                }
            }
        }
        clock_t end_time = clock();
        timeCostEachAssemblyBudget[r] = end_time-start_time+1;
    }
    
    cout<<"finalizing"<<endl;
    for(int r=1; r<=k; r++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
            for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
            // for(int submask = 1; submask < mask; submask++){
                for(int subr=1; subr<=r; subr++){
                    if(0 == r-subr) continue;
                    if(dp[mask][sink][r] > dp[submask][sink][subr]+dp[mask^submask][sink][r-subr]){
                        dp[mask][sink][r] = dp[submask][sink][subr]+dp[mask^submask][sink][r-subr];    
#ifdef RECOVER_MSPATREE
                        triTree[mask][sink][r] = prevTriState(sink, subr, submask, r-subr, mask^submask);
#endif
                    }
                }
            }
            for(int u=1; u<v+1; u++){
                if(dp[mask][sink][r] > dp[mask][u][r]+distanceQuery(u, sink)){
                    dp[mask][sink][r] = dp[mask][u][r]+distanceQuery(u, sink);
#ifdef RECOVER_MSPATREE
                    triTree[mask][sink][r] = prevTriState(u, r, mask, 0, 0);
#endif
                }
            }
        }
        clock_t end_time = clock();
        timeCostEachAssemblyBudget[r] += end_time-start_time;
    }
    /* cout<<" sink to roots "<<endl; */
    /* cout<<" ++++++++"<<endl; */
    /* cout<<" print the result panel "<<endl; */
    /* for(int k = (1<<sNum)-1; k>=0; k--){ */
    /*     cout<<" mask "<<k<<" maskCnt "<<bitCnt(k)<<" bitset "<<bitset<5> (k)<<endl; */
    /*     for(int i = 1; i<v+1; i++){ */
    /*         cout<<" nodeid "<<i; */
    /*         for(int j = 0; j<=3; j++){ */
    /*             cout<<" "<<dp[k][i][j]<<" "; */
    /*         } */
    /*         cout<<endl; */
    /*     } */
    /*     cout<<endl; */
    /* } */
    /* cout<<endl; */

    // for(int i=0;i<(1<<sNum);i++){
    //     cout<<i<<endl;
    // }
    // for(int r=0; r<k+1; r++){
    //     cout<<"r="<<r<<endl;
    //     for(int i=0; i<(1<<sNum); i++){
    //         cout<<"mask="<<i<<endl;
    //         for(int j=0; j<v+1; j++){
                
    //             Quad ele = preState[i][j][r];
    //             cout<<"["<<ele.nodeid<<", "<<ele.s1<<", "<<ele.r1<<", "<<ele.s2<<", "<<ele.r2<<"] "<<" ";
    //             if(j%5==0){
    //                 cout<<endl;
    //             }
    //         }
    //         cout<<endl;
    //     }
    // }

    // fstream fstm_comp_adv;
    // fstm_comp_adv.open("result/comp_bsc.res", ios::out); 
    // for(int mask=1; mask< 1<<sNum; mask++){
    //     for(int i=1; i<v+1; i++){
    //         if(i+1 == v+1) continue;
    //         fstm_comp_adv<<dp[mask][i][1]<<" ";
    //     }
    //     fstm_comp_adv<<endl;
    // }
    // fstm_comp_adv.close();

    for(int i=0; i<timeCostEachAssemblyBudget.size(); i++){
        if(timeCostEachAssemblyBudget[i]!=0){
            cout<<endl<<"budget: "<< i<<endl;
            cout<<"time cost: "<<(double)timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC<<endl;
            cout<<"minDist:  "<<dp[full_mask][sink][i]<<endl;
            res50Dist[i].push_back(dp[full_mask][sink][i]);
            res50Time[i].push_back((double)(timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC));
        }
    }
    
#ifdef RECOVER_MSPATREE
    dfsRTree_BasicDP(triTree, sink, k, full_mask);
#endif 
    for (int i = 0; i < bmsk; i++) {
        for (int j = 0; j < vnum; j++)
        {
            free(dp[i][j]);
        }
        free(dp[i]);
    }
    free(dp);
#ifdef RECOVER_MSPATREE
    for (int i = 0; i < bmsk; i++) {
        for (int j = 0; j < vnum; j++)
        {
            delete[] triTree[i][j];
        }
        delete[] triTree[i];
    }
    delete[] triTree;
#endif
}

//terminals={srcs}∪{sink}, k is the budget of assembly points
//Warning of INF. Given the fact that we have the prior knowledge of max value for the limited parameters of advDP,
//we init the max value of three quater of the MAX_INT. By doing so, we can avoid data type over flow on int type.
void algo6_AdvDP(vector<int>& terminals, int k){
    //minDist shoud be INF, so that it can decrease in iterations.
    int minDist = INF;
    //num of terminals
    int t = terminals.size();
    //num of src
    int sNum = t-1;
    //full_mask
    int full_mask = (1<<sNum)-1;
    //nums of nodes;
    int v = GraphSize;
    vector<int> Srcs(terminals.begin(), terminals.end()-1);
    int sink = terminals.back();
    cout<<"Src: ";
    for(auto ele:Srcs){
        cout<<ele<<" ";
    }
    cout<<"sink: "<<sink<<endl;
    priority_queue<pair<int, int> > pq;

    //[mask][nodid][rpCnt][ifRP]
    //mask 1~2^t-1 ; nodeid 1~n, nodeid starts from 1; rpCnt 0~r; ifRP 0 1;
    // vector<vector<vector<vector<int> > > > dpQ(1<<sNum, \
    //              vector<vector<vector<int> > >(v+1, \
    //                      vector<vector<int> > (bitCnt((1<<sNum)-1)+1, \
    //                                vector<int>(2, INF))));
    // we manually apply for heap memory to reduce cost.
    int bmsk = 1<<sNum;
    int vnum = v+1;
    int bgt = bitCnt((1<<sNum)-1)+1;
    int two = 2;
    int bmsk = 1 << sNum;
    // int ****dpQ = new int***[bmsk];
    // for (int i = 0; i < bmsk; i++) {
    //     dpQ[i] = new int**[vnum];
    //     for (int j = 0; j < vnum; j++) {
    //         dpQ[i][j] = new int*[bgt];
    //         for (int kth = 0; kth < bgt; kth++) {
    //             dpQ[i][j][kth] = new int[two];
    //         }
    //     }
    // }

    // // Initialize the array
    // for (int i = 0; i < (1 << sNum); i++) {
    //     for (int j = 0; j < v + 1; j++) {
    //         for (int kth = 0; kth < bgt; kth++) {
    //             for (int l = 0; l < 2; l++) {
    //                 dpQ[i][j][kth][l] = TRIQUATINF; // Assuming TRIQUATINF is a defined constant
    //             }
    //         }
    //     }
    // }

    // // Don't forget to free the memory when done
    // for (int i = 0; i < bmsk; i++) {
    //     for (int j = 0; j < vnum; j++) {
    //         for (int kth = 0; kth < bgt; kth++) {
    //             delete[] dpQ[i][j][kth]; // Free the innermost array
    //         }
    //         delete[] dpQ[i][j]; // Free the 3rd level array
    //     }
    //     delete[] dpQ[i]; // Free the 2nd level array
    // }
    // delete[] dpQ; // Free the top-level array
    int ****dpQ = (int ****)malloc(bmsk * sizeof(int ***));
    for(int i=0; i<bmsk; i++){
        dpQ[i] = (int ***)malloc(vnum * sizeof(int **));
        for(int j=0; j<vnum; j++){
            dpQ[i][j] = (int **)malloc(bgt * sizeof(int *));
            for(int kth=0; kth<bgt; kth++){
                dpQ[i][j][kth] = (int *)malloc(two * sizeof(int));
            }
        }
    }
    for(int i=0; i<(1<<sNum); i++){
        for(int j=0; j<v+1; j++){
            for(int kth=0; kth<(bgt); kth++){
                for(int l=0; l<2; l++){
                    dpQ[i][j][kth][l] = TRIQUATINF;
                }
            }
        }
    }
    // vector<int> vis(v+1, 0);
    for(int i=0; i<Srcs.size(); i++){
        dpQ[1<<i][Srcs[i]][0][0] = 0;
        dpQ[1<<i][Srcs[i]][1][1] = 0;
    }
#ifdef RECOVER_MSPATREE
    //CAUTION: It takes more than twice of the memory if retrieving MSPA tree turns on.
    //state init
    prevState**** rTree = new prevState***[bmsk];
    for(int i=0; i<bmsk; i++){
        rTree[i] = new prevState**[vnum];
        for(int j=0; j<vnum; j++){
            rTree[i][j] = new prevState*[bgt];
            for(int kth=0; kth<bgt; kth++){
                rTree[i][j][kth] = new prevState[two];
            }
        }
    }
#endif    
    // Dijkstra to init distance from src to the assembly point when k=1,
    // because src to its nearest assembly(include sink) is shortest path and therefore equal to k=0.
    // We do not exclude dist to sink here, 
    // because if srcs go direct to the sink is the cheapest way, 
    // the only assembly point is inexplicitly included on any shortest path from srcs to sink.
    // Later other states will be computed in the next part.
    for(int i=0; i<Srcs.size(); i++){
        vector<int> vis(v+1, 0);
        priority_queue<pair<int, int> > pq_local;
        pq_local.emplace(make_pair(0, Srcs[i]));
        while(!pq_local.empty()){
            pair<int, int> pii = pq_local.top();
            int u = pii.second;
            pq_local.pop();
            if(vis[u]) continue;
            vis[u] = 1;
            for(int kth=0; kth<Degree[u]; kth++){
                int nei = Neighbor[u][kth];
                if(dpQ[1<<i][nei][0][0] > dpQ[1<<i][u][0][0]+Weight[u][kth]){
                    dpQ[1<<i][nei][0][0] = dpQ[1<<i][u][0][0]+Weight[u][kth];
                    dpQ[1<<i][nei][1][1] = dpQ[1<<i][nei][0][0];
#ifdef RECOVER_MSPATREE
                    //init states have no prevStates.
                    rTree[1<<i][nei][1][1] = prevState(0, 0, 0, 0, 0, 0, 0);
                    rTree[1<<i][nei][0][0] = prevState(0, 0, 0, 0, 0, 0, 0);
#endif
                    pq_local.emplace(make_pair(-dpQ[1<<i][nei][1][1], nei));
                }
            }
        }
    } 

    cout<<"init fin"<<endl;
    //tree merge & tree grow  
    vector<clock_t> timeCostEachAssemblyBudget(k+2, 0);
    for(int rcnt=1; rcnt<=k; rcnt++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
            //@caution nodeid starts form 1
            for(int nodeid = 1; nodeid <v+1; nodeid++){
                if(nodeid == sink) continue;
                //tree merge
                //enuerate all subset of mask. e.g. (0010) has no subset.
                for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
                // for(int submask=1; submask<mask; submask++){
                    for(int subr=1; subr<=rcnt; subr++){
                        if(dpQ[mask][nodeid][rcnt][1] >  dpQ[submask][nodeid][subr][1]+dpQ[submask^mask][nodeid][(rcnt+1)-subr][1]){
                            dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], \
                                                             dpQ[submask][nodeid][subr][1]+dpQ[submask^mask][nodeid][(rcnt+1)-subr][1]);
#ifdef RECOVER_MSPATREE
                            //tree merge operation merges at the root by two different states at the same root.
                            rTree[mask][nodeid][rcnt][1] = prevState(nodeid, subr, submask, 1, (rcnt+1)-subr, mask^submask, 1);
#endif
                        }
                    }
                }
                //tree grow
                for(int nei=0; nei<Degree[nodeid]; nei++){
                    int u = Neighbor[nodeid][nei];
                    //sink cannot be marked as a assembly as it will never consume budget
                    if(u == sink) continue;
                    // already  initialized
                    if((rcnt-1) <1 ) continue;
                    if(dpQ[mask][nodeid][rcnt][1] >  dpQ[mask][u][rcnt-1][0]+Weight[nodeid][nei]){
                        /* if(rcnt-1 == 0) continue; */
                        dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], \
                                                         dpQ[mask][u][rcnt-1][0]+Weight[nodeid][nei]);
#ifdef RECOVER_MSPATREE
                         //the tree grow keeps one preve states. 0 -> 1
                        rTree[mask][nodeid][rcnt][1] = prevState(u, rcnt-1, mask, 0, 0, 0, 0);
#endif
                    }
                    if(dpQ[mask][nodeid][rcnt][1] >  dpQ[mask][u][rcnt-1][1]+Weight[nodeid][nei]){ 
                        dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1],  
                                                         dpQ[mask][u][rcnt-1][1]+Weight[nodeid][nei]); 
#ifdef RECOVER_MSPATREE
                        //1->1
                        rTree[mask][nodeid][rcnt][1] = prevState(u, rcnt-1, mask, 1, 0, 0, 0);
#endif
                    } 
                }
                pq.emplace(make_pair(-dpQ[mask][nodeid][rcnt][1], nodeid));
                while(!pq.empty()){
                    pair<int, int> vk = pq.top();
                    int v = vk.second;
                    pq.pop();
                    for(int nei=0; nei<Degree[v]; nei++){
                        int u = Neighbor[v][nei];
                        if(u == sink) continue;
                        if(dpQ[mask][u][rcnt][0] > min(dpQ[mask][v][rcnt][0]+Weight[v][nei], dpQ[mask][v][rcnt][1]+Weight[v][nei])){
                            dpQ[mask][u][rcnt][0] = min(dpQ[mask][u][rcnt][0], \
                                                        min(dpQ[mask][v][rcnt][0]+Weight[v][nei], \
                                                            dpQ[mask][v][rcnt][1]+Weight[v][nei]));
#ifdef RECOVER_MSPATREE
                            // 0->0 or 1->0
                            if(dpQ[mask][v][rcnt][0]+Weight[v][nei] >  dpQ[mask][v][rcnt][1]+Weight[v][nei]){
                               rTree[mask][u][rcnt][0] = prevState(v, rcnt, mask, 1, 0, 0, 0);
                            }else{ 
                               rTree[mask][u][rcnt][0] = prevState(v, rcnt, mask, 0, 0, 0, 0);
                            } 
#endif
                            pq.emplace(make_pair(-dpQ[mask][u][rcnt][0], u));
                        }
                    }
                }
                
            }
        }
        clock_t end_time= clock();
        // add 1e-6 to the time record for test on example.gr.
        timeCostEachAssemblyBudget[rcnt] = end_time-start_time+1;
    }
    cout<<"finallizing"<<endl;
    // destination t does not consume any budget, therefore rcnt can be as large as k
    for(int rcnt=1; rcnt<=k; rcnt++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
                for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
                // for(int submask = 1; submask < mask; submask++){
                    for(int subr=1; subr<=rcnt; subr++){
                         if(dpQ[mask][sink][rcnt][1] >  dpQ[submask][sink][subr][1]+dpQ[submask^mask][sink][rcnt-subr][1]){ 
                             dpQ[mask][sink][rcnt][1] = min(dpQ[mask][sink][rcnt][1], 
                                                            dpQ[submask][sink][subr][1]+dpQ[submask^mask][sink][rcnt-subr][1]); 

#ifdef RECOVER_MSPATREE
                            //The states grow to sink can merge at the root.
                            rTree[mask][sink][rcnt][1] = prevState(sink, subr, submask, 1, rcnt-subr, mask^submask, 1);
#endif
                         } 
                    }
                }
                for(int nei=0; nei<Degree[sink]; nei++){
                    int vtosink = Neighbor[sink][nei];
                    if( dpQ[mask][sink][rcnt][1] > min(dpQ[mask][vtosink][rcnt][0]+Weight[sink][nei], dpQ[mask][vtosink][rcnt][1]+Weight[sink][nei])){
                        dpQ[mask][sink][rcnt][1] = min(dpQ[mask][sink][rcnt][1], 
                                                       min(dpQ[mask][vtosink][rcnt][1]+Weight[sink][nei], 
                                                           dpQ[mask][vtosink][rcnt][0]+Weight[sink][nei]));
#ifdef RECOVER_MSPATREE
                       //0->1 or 1->1
                        if(dpQ[mask][vtosink][rcnt][1]+Weight[sink][nei] > dpQ[mask][vtosink][rcnt][0]+Weight[sink][nei]){ 
                           rTree[mask][sink][rcnt][1] = prevState(vtosink, rcnt, mask, 0, 0, 0, 0);
                        }else{ 
                           rTree[mask][sink][rcnt][1] = prevState(vtosink, rcnt, mask, 1, 0, 0, 0);
                        } 
#endif
                    }
                }
        }
        clock_t end_time = clock();
        timeCostEachAssemblyBudget[rcnt] += end_time-start_time;
    }
    // for(int mask=1; mask< 1<<sNum; mask++){
    //     for(int i=1; i<v+1; i++){
    //         cout<<dpQ[mask][i][1][1]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"-------"<<endl;
    // for(int rcnt=1; rcnt<=k; rcnt++){
    //     for(int i=1; i<v+1; i++){
    //         cout<<dpQ[full_mask][sink][rcnt][1]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // fstream fstm_comp_adv;
    // fstm_comp_adv.open("result/comp_adv.res", ios::out); 
    // for(int mask=1; mask< 1<<sNum; mask++){
    //     for(int i=1; i<v+1; i++){
    //         if(i+1 == v+1) continue;
    //         fstm_comp_adv<<dpQ[mask][i][1][1]<<" ";
    //     }
    //     fstm_comp_adv<<endl;
    // }
    // fstm_comp_adv.close();
    //res to file
    for(int i=0; i<timeCostEachAssemblyBudget.size(); i++){
        if(timeCostEachAssemblyBudget[i]!=0){
            cout<<endl<<"budget: "<< i<<endl;
            cout<<"time cost: "<<(double)timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC<<endl;
            cout<<"minDist:  "<<dpQ[full_mask][sink][i][1]<<endl;
            res50Dist[i].push_back(dpQ[full_mask][sink][i][1]);
            res50Time[i].push_back((double)(timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC));
        }
    }
#ifdef RECOVER_MSPATREE
    //retrieve assembly points
    dfsRetrieveTree(rTree, sink, k, full_mask, 1);
#endif
    //do not forget to free!!!
    for (int i = 0; i < bmsk; i++) {
        for (int j = 0; j < vnum; j++)
        {
            for (int kth = 0; kth < bgt; kth++)
            {
                free(dpQ[i][j][kth]);
            }
            free(dpQ[i][j]);
        }
        free(dpQ[i]);
    }
    free(dpQ);
#ifdef RECOVER_MSPATREE
    //do not forget to free!!!
    for (int i = 0; i < bmsk; i++) {
        for (int j = 0; j < vnum; j++)
        {
            for (int kth = 0; kth < bgt; kth++)
            {
                delete[] rTree[i][j][kth];
            }
            delete[] rTree[i][j];
        }
        delete[] rTree[i];
    }
    delete[] rTree;
#endif
}
//@Note: candidates do not include the sink vertex.
void algo6_AdvDP_given_candidates(unordered_set<int>& candidates, vector<int>& terminals, int k ){
    //minDist shoud be INF, so that it can decrease in iterations.
    int minDist = INF;
    //num of terminals
    int t = terminals.size();
    //num of src
    int sNum = t-1;
    //full_mask
    int full_mask = (1<<sNum)-1;
    //nums of nodes;
    int v = GraphSize;
    //num of candidates
    int candNum = candidates.size();
    vector<int> Srcs(terminals.begin(), terminals.end()-1);
    int sink = terminals.back();
    cout<<"src: ";
    for(auto ele:Srcs){
        cout<<ele<<" ";
    }
    cout<<endl;
    cout<<"sink: "<<sink<<endl;
    priority_queue<pair<int, int> > pq;

    //[mask][nodid][rpCnt][ifRP]
    int bmsk = 1<<sNum;
    int vnum = v+1;
    int bgt = bitCnt((1<<sNum)-1)+1;
    int two = 2;
    int ****dpQ = (int ****)malloc(bmsk * sizeof(int ***));
    for(int i=0; i<bmsk; i++){
        dpQ[i] = (int ***)malloc(vnum * sizeof(int **));
        for(int j=0; j<vnum; j++){
            dpQ[i][j] = (int **)malloc(bgt * sizeof(int *));
            for(int kth=0; kth<bgt; kth++){
                dpQ[i][j][kth] = (int *)malloc(two * sizeof(int));
            }
        }
    }
    // int dpQ[1<<sNum][v+1][bitCnt((1<<sNum)-1)+1][2]; 
    for(int i=0; i<(1<<sNum); i++){
        for(int j=0; j<v+1; j++){
            for(int kth=0; kth<(bitCnt((1<<sNum)-1)+1); kth++){
                for(int l=0; l<2; l++){
                    dpQ[i][j][kth][l] = INF;
                }
            }
        }
    }
    vector<int> vis(v+1, 0);
    for(int i=0; i<Srcs.size(); i++){
        dpQ[1<<i][Srcs[i]][0][0] = 0;
        dpQ[1<<i][Srcs[i]][1][1] = 0;
    }
    // The init phase stays the same.
    // The only difference lies in the limitation on non-candidate vertices. 
    // As only the candidates are legal as assembly points, other vertices can only conduct tree grow operations.
    for(int i=0; i<Srcs.size(); i++){
        vector<int> vis(v+1, 0);
        priority_queue<pair<int, int> > pq;
        pq.emplace(make_pair(0, Srcs[i]));
        // vis[Srcs[i]] = 1;
        while(!pq.empty()){
            pair<int, int> pii = pq.top();
            int u = pii.second;
            pq.pop();
            if(vis[u]) continue;
            vis[u] = 1;
            for(int kth=0; kth<Degree[u]; kth++){
                int nei = Neighbor[u][kth];
                // if (vis[nei]) continue;
                // vis[nei] = 1;
                if(dpQ[1<<i][nei][0][0] > dpQ[1<<i][u][0][0]+Weight[u][kth]){
                    dpQ[1<<i][nei][0][0] = dpQ[1<<i][u][0][0]+Weight[u][kth];
                    dpQ[1<<i][nei][1][1] = dpQ[1<<i][nei][0][0];
                    pq.emplace(make_pair(-dpQ[1<<i][nei][1][1], nei));
                }
            }
        }
    } 
    cout<<"init fin"<<endl;
    //tree merge & tree grow  
    //Non-candate vertices are only allowed for tree grow operations.
    vector<clock_t> timeCostEachAssemblyBudget(k+2, 0);
    for(int rcnt=1; rcnt<=k; rcnt++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
            //@caution nodeid starts from 1.
            for(int nodeid = 1; nodeid < v+1; nodeid++){
                if(nodeid == sink) continue;
                // tree merge
                if(candidates.find(nodeid) != candidates.end()){
                    // If not a member of candidate set, then tree merge operations will not conduct on this vertex.
                    for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
                    // for (int submask = 1; submask < mask; submask++){
                        for (int subr = 1; subr <= rcnt; subr++){
                            if (dpQ[mask][nodeid][rcnt][1] > dpQ[submask][nodeid][subr][1] + dpQ[submask ^ mask][nodeid][(rcnt + 1) - subr][1]){
                                dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1],
                                                                 dpQ[submask][nodeid][subr][1] + dpQ[submask ^ mask][nodeid][(rcnt + 1) - subr][1]);
                            }
                        }
                    }
                    // tree grow
                    // If a vertex is an assembly point, then it must be a member of candidate set.
                    for (int nei = 0; nei < Degree[nodeid]; nei++){
                        int u = Neighbor[nodeid][nei];
                        // sink cannot be marked as a assembly as it will never consume budget
                        if (u == sink) continue;
                        // already  initialized
                        if ((rcnt - 1) < 1) continue;
                        if (dpQ[mask][nodeid][rcnt][1] > dpQ[mask][u][rcnt - 1][0] + Weight[nodeid][nei]){
                            /* if(rcnt-1 == 0) continue; */
                            dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], \
                                                             dpQ[mask][u][rcnt - 1][0] + Weight[nodeid][nei]);
                        }
                        if (dpQ[mask][nodeid][rcnt][1] > dpQ[mask][u][rcnt - 1][1] + Weight[nodeid][nei]){
                            dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], \
                                                             dpQ[mask][u][rcnt - 1][1] + Weight[nodeid][nei]);
                        }
                    }
                }
                pq.emplace(make_pair(-dpQ[mask][nodeid][rcnt][1], nodeid));
                // cout<<"mask: "<<mask<<" nodeid: "<<nodeid<<" rcnt: "<<rcnt<<" val: "<<dpQ[mask][nodeid][rcnt][1]<<endl;
                while(!pq.empty()){
                    pair<int, int> vk = pq.top();
                    int v = vk.second;
                    pq.pop();
                    for(int nei=0; nei<Degree[v]; nei++){
                        int u = Neighbor[v][nei];
                        if(u == sink) continue;
                        if(dpQ[mask][u][rcnt][0] > min(dpQ[mask][v][rcnt][0]+Weight[v][nei], dpQ[mask][v][rcnt][1]+Weight[v][nei])){
                            dpQ[mask][u][rcnt][0] = min(dpQ[mask][u][rcnt][0], \
                                                        min(dpQ[mask][v][rcnt][0]+Weight[v][nei], \
                                                            dpQ[mask][v][rcnt][1]+Weight[v][nei]));
                            pq.emplace(make_pair(-dpQ[mask][u][rcnt][0], u));
                        }
                    }
                }
                
            }
        }
        clock_t end_time= clock();
        timeCostEachAssemblyBudget[rcnt] = end_time-start_time+1;
        // cout<<"#AP: "<<rcnt<<" time: "<<timeCostEachAssemblyBudget[rcnt]<<endl;
    }
    cout<<"finallizing"<<endl;
    // destination t does not consume any budget, therefore rcnt can be as large as k
    for(int rcnt=1; rcnt<=k; rcnt++){
        clock_t start_time = clock();
        for(int mask = 1; mask< 1<<sNum; mask++){
                for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
                // for(int submask = 1; submask < mask; submask++){
                    for(int subr=1; subr<=rcnt; subr++){
                         if(dpQ[mask][sink][rcnt][1] >  dpQ[submask][sink][subr][1]+dpQ[submask^mask][sink][rcnt-subr][1]){ 
                             dpQ[mask][sink][rcnt][1] = min(dpQ[mask][sink][rcnt][1], 
                                                            dpQ[submask][sink][subr][1]+dpQ[submask^mask][sink][rcnt-subr][1]); 

                         } 
                    }
                }
                for(int nei=0; nei<Degree[sink]; nei++){
                    int vtosink = Neighbor[sink][nei];
                    if( dpQ[mask][sink][rcnt][1] > min(dpQ[mask][vtosink][rcnt][0]+Weight[sink][nei], dpQ[mask][vtosink][rcnt][1]+Weight[sink][nei])){
                        dpQ[mask][sink][rcnt][1] = min(dpQ[mask][sink][rcnt][1], 
                                                       min(dpQ[mask][vtosink][rcnt][1]+Weight[sink][nei], 
                                                           dpQ[mask][vtosink][rcnt][0]+Weight[sink][nei]));
                    }
                }
        }
        clock_t end_time = clock();
        timeCostEachAssemblyBudget[rcnt] += end_time-start_time;
        // cout<<"#AP: "<<rcnt<<" time: "<<timeCostEachAssemblyBudget[rcnt]<<endl;
    }
    // for(int mask=1; mask< 1<<sNum; mask++){
    //     for(int i=1; i<v+1; i++){
    //         cout<<dpQ[mask][i][1][1]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    //res to file
    for(int i=0; i<timeCostEachAssemblyBudget.size(); i++){
        if(timeCostEachAssemblyBudget[i]!=0){
            cout<<endl<<"budget: "<< i<<endl;
            cout<<"time cost: "<<(double)timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC<<endl;
            cout<<"minDist:  "<<dpQ[full_mask][sink][i][1]<<endl;
            res50Dist[i].push_back(dpQ[full_mask][sink][i][1]);
            res50Time[i].push_back((double)(timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC));
        }
    }
    //do not forget to free!!!
    for (int i = 0; i < bmsk; i++) {
        for (int j = 0; j < vnum; j++)
        {
            for (int kth = 0; kth < bgt; kth++)
            {
                free(dpQ[i][j][kth]);
            }
            free(dpQ[i][j]);
        }
        free(dpQ[i]);
    }
    free(dpQ);
}
// //terminals={srcs}+{sink}, k is the budget of assembly points
// void algo6_AdvDP_Social(vector<int>& terminals, int k){
//     //minDist shoud be INF, so that it can decrease in iterations.
//     int minDist = INF;
//     //num of terminals
//     int t = terminals.size();
//     //num of src
//     int sNum = t-1;
//     //full_mask
//     int full_mask = (1<<sNum)-1;
//     //nums of nodes;
//     int v = GraphSize;
//     vector<int> Srcs(terminals.begin(), terminals.end()-1);
//     int sink = terminals.back();
//     cout<<"Src "<<endl;
//     for(auto ele:Srcs){
//         cout<<ele<<" ";
//     }
//     cout<<" sink "<<sink<<endl;
//     priority_queue<pair<int, int> > pq;

//     //[mask][nodid][rpCnt][ifRP]
//     //mask 1~2^t-1 ; nodeid 1~n, nodeid starts from 1; rpCnt 0~r; ifRP 0 1;
//     vector<vector<vector<vector<int> > > > dpQ(1<<sNum, \
//                  vector<vector<vector<int> > >(v+1, \
//                          vector<vector<int> > (bitCnt((1<<sNum)-1)+1, \
//                                    vector<int>(2, INF))));
//     vector<int> vis(v+1, 0);
//     for(int i=0;i<Srcs.size(); i++){
//         for(int j=1; j<v+1; j++){
//            dpQ[1<<i][j][0][0] = pll.QueryDistance(Srcs[i], j); 
//         }
//     }
//     for(int i=0; i<Srcs.size(); i++){
//         dpQ[1<<i][Srcs[i]][0][0] = 0;
//         dpQ[1<<i][Srcs[i]][1][1] = 0;
//     }
// #ifdef RECOVER_MSPATREE
//     Pent pt;
//      for(int i=0;i<(1<<sNum);i++){ 
//          pre.push_back(vector<vector<vector<Pent> >> (v+1, vector<vector<Pent>>(bitCnt((1<<sNum)-1)+1, vector<Pent>(2, pt)))); 
//      } 
//      for(int i=0;i<(1<<sNum); i++){ 
//          foot.push_back(vector<vector<vector<int> >> (v+1, vector<vector<int>>(bitCnt((1<<sNum)-1)+1, vector<int>(2, 0)))); 
//      } 
// #endif    
//     for(int i=0; i<Srcs.size(); i++){
//         vector<int> vis(v+1, 0);
//         priority_queue<pair<int, int> > pq;
//         pq.emplace(make_pair(0, Srcs[i]));
//         while(!pq.empty()){
//             pair<int, int> pii = pq.top();
//             int u = pii.second;
//             pq.pop();
//             if(vis[u]) continue;
//             vis[u] = 1;
//             for(int k=0; k<Degree[u]; k++){
//                 int nei = Neighbor[u][i];
//                 if(dpQ[1<<i][nei][1][1] > dpQ[1<<i][u][0][0]+Weight[u][k]){
//                     dpQ[1<<i][nei][1][1] = dpQ[1<<i][u][0][0]+Weight[u][k];
//                     //recover MSPA tree
// #ifdef RECOVER_MSPATREE
//                      pre[1<<i][nei][1][1] = Pent(u, (1<<i), 1, (1<<i), 1, 0);
// #endif
//                     pq.emplace(make_pair(-dpQ[1<<i][nei][1][1], nei));
//                 }
//             }
//         }
//     } 
//     cout<<"init fin"<<endl;
//     //tree merge & tree grow  
//     vector<clock_t> timeCostEachAssemblyBudget(k+2, 0);
//     for(int rcnt=1; rcnt<=k; rcnt++){
//         clock_t start_time = clock();
//         for(int mask = 1; mask< 1<<sNum; mask++){
//             //@caution nodeid starts form 1
//             for(int nodeid = 1; nodeid <v+1; nodeid++){
//                 if(nodeid == sink) continue;
//                 for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
//                     for(int subr=1; subr<=rcnt; subr++){
//                         if(dpQ[mask][nodeid][rcnt][1] >  dpQ[submask][nodeid][subr][1]+dpQ[submask^mask][nodeid][(rcnt+1)-subr][1]){
//                             dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], \
//                                                              dpQ[submask][nodeid][subr][1]+dpQ[submask^mask][nodeid][(rcnt+1)-subr][1]);
// #ifdef RECOVER_MSPATREE
//                              pre[mask][nodeid][rcnt][1] = Pent(nodeid, submask, subr, mask^submask, (rcnt+1)-subr, 1); 
// #endif
//                         }
//                     }
//                 }
//                 for(int nei=0; nei<Degree[nodeid]; nei++){
//                     int u = Neighbor[nodeid][nei];
//                     if(dpQ[mask][nodeid][rcnt][1] >  dpQ[mask][u][rcnt-1][0]+Weight[nodeid][nei]){
//                         /* if(rcnt-1 == 0) continue; */
//                         dpQ[mask][nodeid][rcnt][1] = min(dpQ[mask][nodeid][rcnt][1], dpQ[mask][u][rcnt-1][0]+Weight[nodeid][nei]);
//                         //recover MSPA tree
// #ifdef RECOVER_MSPATREE
//                      pre[mask][nodeid][rcnt][1] = Pent(u, mask, rcnt-1, mask, rcnt-1, 0);
// #endif
//                     }
//                     //sink cannot be marked as a assembly as it will never consume budget
//                     if(u == sink) continue;
//                     //buget should be greater than 1 if it is a assembly point
//                     if((rcnt-1) <1 ) continue;
//                 }
//                 pq.emplace(make_pair(-dpQ[mask][nodeid][rcnt][1], nodeid));
//                 while(!pq.empty()){
//                     pair<int, int> vk = pq.top();
//                     int v = vk.second;
//                     pq.pop();
//                     for(int nei=0; nei<Degree[v]; nei++){
//                         int u = Neighbor[v][nei];
//                         if(dpQ[mask][u][rcnt][0] > min(dpQ[mask][v][rcnt][0]+Weight[v][nei], dpQ[mask][v][rcnt][1]+Weight[v][nei])){
//                             dpQ[mask][u][rcnt][0] = min(dpQ[mask][u][rcnt][0], \
//                                                         min(dpQ[mask][v][rcnt][0]+Weight[v][nei], \
//                                                             dpQ[mask][v][rcnt][1]+Weight[v][nei]));
// #ifdef RECOVER_MSPATREE
//                              if(dpQ[mask][v][rcnt][0]+Weight[v][nei] >  dpQ[mask][v][rcnt][1]+Weight[v][nei]){
//                                  pre[mask][u][rcnt][0] = Pent(v, mask, rcnt, mask, rcnt, 1); 
//                              }else{ 
//                                  pre[mask][u][rcnt][0] = Pent(v, mask, rcnt, mask, rcnt, 0); 
//                              } 
// #endif
//                             pq.emplace(make_pair(-dpQ[mask][u][rcnt][0], u));
//                         }
//                     }
//                 }
                
//             }
//         }
//         clock_t end_time= clock();
//         timeCostEachAssemblyBudget[rcnt] = end_time-start_time;
//     }
//     cout<<"finallizing"<<endl;
//     for(int rcnt=1; rcnt<=k; rcnt++){
//         clock_t start_time = clock();
//         for(int mask = 1; mask< 1<<sNum; mask++){
//                 for(int submask = (mask-1)&mask; submask>0; submask = (submask-1)&mask){
//                     for(int subr=1; subr<=rcnt; subr++){
//                          if(dpQ[mask][sink][rcnt][1] >  dpQ[submask][sink][subr][0]+dpQ[submask^mask][sink][rcnt-subr][0]){ 
//                              dpQ[mask][sink][rcnt][0] = min(dpQ[mask][sink][rcnt][0], 
//                                                             dpQ[submask][sink][subr][0]+dpQ[submask^mask][sink][rcnt-subr][0]); 

// #ifdef RECOVER_MSPATREE
//                              pre[mask][sink][rcnt][0] = Pent(sink, submask, subr, mask^submask, rcnt-subr, 0); 
// #endif
//                          } 
//                     }
//                 }
//                 for(int i=1; i<v+1; i++){
//                     if( dpQ[mask][sink][rcnt][0] > min(dpQ[mask][i][rcnt][1]+pll.QueryDistance(i, sink), dpQ[mask][i][rcnt][0]+pll.QueryDistance(i, sink))){
//                         dpQ[mask][sink][rcnt][0] = min(dpQ[mask][sink][rcnt][0], 
//                                                        min(dpQ[mask][i][rcnt][1]+pll.QueryDistance(i, sink), 
//                                                            dpQ[mask][i][rcnt][0]+pll.QueryDistance(i, sink)));
// #ifdef RECOVER_MSPATREE
//                          if(dpQ[mask][i][rcnt][1]+pll.QueryDistance(i, sink) > dpQ[mask][i][rcnt][0]+pll.QueryDistance(i, sink)){ 
//                              pre[mask][sink][rcnt][0] = Pent(i, mask, rcnt, mask, rcnt, 0); 
//                          }else{ 
//                              pre[mask][sink][rcnt][0] = Pent(i, mask, rcnt, mask, rcnt, 1); 
//                          } 
// #endif
//                     }
//                 }
//         }
//         clock_t end_time = clock();
//         timeCostEachAssemblyBudget[rcnt] += end_time-start_time;
//     }
//     //res to file
//     for(int i=0; i<timeCostEachAssemblyBudget.size(); i++){
//         if(timeCostEachAssemblyBudget[i]!=0){
//             cout<<"number of budeget i "<< i<<" time cost "<<(double)timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC<<endl;
//             cout<<"minDist  "<<dpQ[full_mask][sink][i][0]<<endl;
//             res50Dist[i].push_back(dpQ[full_mask][sink][i][0]);
//             res50Time[i].push_back((double)(timeCostEachAssemblyBudget[i]/CLOCKS_PER_SEC));
//         }
//     }
// }


long long algo_simpleGreedy(vector<int>& terminal, int budget, vector<int>& aps,  long long& upDist){

    vector<int> src(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<vector<int>> vvi_times(terminal.size(), vector<int>(GraphSize+1, 0));
    vector<int> totDist(GraphSize+1, 0);
    vector<int> longest_edge_of_src(terminal.size()+1, 0);
    int cliqueSize = terminal.size()+budget;
    vector<vector<pair<int, int>>> newG(cliqueSize, vector<pair<int, int>>(cliqueSize, make_pair(0, 0)));
    queue<int> qi;
    if(budget == 1){
        // queue<int> qi;
        int minDist = 0;
        // tot dist should be no larger than direct to sink
        for (auto ele : src)
        {
            minDist += distanceQuery(ele, sink);
        }
        cout << "upperbound: " << minDist << endl;
        int first_ap = 0;
        int first_ap_cnt = 0;
        for (int i = 1; i < GraphSize + 1; i++)
        {
            if (i == sink)
                continue;
            int tmp_ap_cnt = 0;
            int curDist = distanceQuery(i, sink);
            for (int j = 0; j < src.size(); j++)
            {
                int src2i = distanceQuery(i, src[j]);
                int src2t = distanceQuery(src[j], sink); 
                if (src2i > src2t)
                {
                    curDist += src2t;
                }
                else
                {
                    curDist += src2i;
                    tmp_ap_cnt++;
                }
            }
            if (minDist > curDist)
            {
                first_ap = i;
                minDist = curDist;
                first_ap_cnt = tmp_ap_cnt;
                // cout << "curDist: " << curDist << endl;
                qi.push(i);
            }
        }
        cout << "minDist 1st k: " << minDist << endl;
        for (int i = 0; i < src.size(); i++)
        {
            longest_edge_of_src[i] = min(distanceQuery(sink, src[i]),\
            distanceQuery(first_ap, src[i]));
        }
        int fir_par = 0;
        for(auto ele:longest_edge_of_src){
            fir_par+=ele;
        }
        cout<<"fir par weight: "<<fir_par<<endl;
        aps.push_back(first_ap);
        cout << "1st ap cnt: " << first_ap_cnt << endl;

        if (budget <= 1 || first_ap_cnt == 0)
        {
            return minDist;
        }
    }
    //k>1
    int qisize = qi.size();
    int cnt = 0;
    int kth = 1;
    int kthMinDist = upDist;
    int kthAP = 0;
    while(!qi.empty()){
        int u = qi.front();
        qi.pop();
        cnt++;
        if(cnt == qisize){
            //no valid next ap, break
            if(kthMinDist >= upDist){
                break;
            }
            // upd upperbound...
            upDist = kthMinDist;
            aps.push_back(kthAP);
            // kth ...
            kth++;
            if(kth >= budget){
                break;
            }
        }
        int excess_cnt = 0;
        for(int j=0; j<src.size(); j++){
            int s2u = distanceQuery(src[j], u);
            if(s2u > longest_edge_of_src[j]){
                excess_cnt++;
            }
        }
        //simple pruning, remove unpromising cands.
        if(excess_cnt < 4){
            //possible cand for next ap
            qi.push(u);
            int s2c=0;
            //assgin to nearest 
            for(int i=0; i<src.size(); i++){
                int src2cand = distanceQuery(u, src[i]);
                //assigned to either u
                if(src2cand < longest_edge_of_src[i]){
                    s2c += src2cand;
                //or stay unchanged.
                }else{
                    s2c += longest_edge_of_src[i];
                }
            }
            cout<<"s2c: "<<s2c<<endl;
            // make clique: src, aps, u
            aps.push_back(u);
            fillG(aps, sink, newG);
            //prim
            int c2t = prim(newG, aps.size()+1);
            aps.pop_back();
            //comp if valid
            if(c2t+s2c < kthMinDist){
                kthMinDist = c2t+s2c;
                kthAP = u;
            }
        }
    }

    return upDist;
}
// loop over all vertices and find next possible ap. 
// O(kn), it is ugly slow.
int algo_greedy(vector<int>& terminal, vector<int>& aps, int budget){
    // init
    vector<int> src(terminal.begin(), terminal.end()-1);
    int cliqueSize = terminal.size()+budget;
    vector<vector<pair<int, int>>> newG(cliqueSize, vector<pair<int, int>>(cliqueSize, make_pair(0, 0)));
    int sink = terminal.back();
    int minDist = 0;
    for(auto ele:src){
        minDist += distanceQuery(ele, sink);
    }
    int ap=0;
    int loopMinDist = minDist;
    while(1){
        int curDist = 0;
        for(int i=1; i<GraphSize+1; i++){
            if(i%10000 == 0){
                cout<<i<<endl;
            }
            if (i == sink)
                continue;
            int curDist = 0; 
            aps.push_back(i);
            for (int j = 0; j < src.size(); j++)
            {
                int src2t = distanceQuery(src[j], sink);
                int src2i = INF; 
                int whichaps = 0;
                for(int kth=0; kth<aps.size(); kth++){
                    int src2loopc =  distanceQuery(aps[kth], src[j]);
                    if(src2i > src2loopc){
                        src2i = src2loopc;
                    }
                }
                if (src2i > src2t)
                {
                    curDist += src2t;
                }
                else
                {
                    curDist += src2i;
                }
            }
            //if 1st part is greater than upbound, no further computation for 2nd part.
            if(curDist > loopMinDist) continue;
            cout<<"src+aps: "<<src.size()+aps.size()<<endl;
            fillG(aps, sink, newG);
            int c2t = prim(newG, src.size()+aps.size());
            curDist += c2t;
            if(loopMinDist > curDist){
                ap = i;
                loopMinDist = curDist;
            }
            aps.pop_back();
        }//end for
        if(minDist > loopMinDist){
            aps.push_back(ap);
            minDist = loopMinDist;
        }else{
            break; 
        }
        cout<<"minDist: "<<minDist<<" kth: "<<aps.size()<<endl;
    }//end while
    return minDist;
}
//worst case O(n^2)
void algo_apx_greedy(vector<int>& terminal, vector<int>& cand, int budget){
    res50Dist.resize(budget+1);
    res50Time.resize(budget+1);
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long ub = 0;
    for(auto ele:src){
        ub+=distanceQuery(ele, sink);
    }
    cout<<"ub: "<<ub<<endl;
    int pos = 0;
    int newStartPos = 0;
    long long minDist = LLINF;
    long long curMinDist = LLINF;
    vector<int> aps;
    //in the search tree, we greedily select the most promising cand each time.
    //Such that we do not need to enum all the possible combinations as in BnB.
    //i starts from 0 because cand starts from zero
    for(int i=0; i<budget; i++){
        auto apx_grdy_start = timerSt();
        for(int j=newStartPos; j<budget; j++){
            aps.push_back(cand[j]);
            vector<vector<pair<int, int>>> newG;
            // auto s2c_st = timerSt();
            unordered_map<int, vector<int>> umpii;
            long long s2c = assign_src2cand(terminal, aps, newG, umpii, 1);
            // auto s2c_nd = timerNd();
            // auto c2t_st = timerSt();
            long long c2t = prim(newG, newG.size());
            // auto c2t_nd = timerNd();
            // double s2c_dur = dur(s2c_nd, s2c_st);
            // double c2t_dur = dur(c2t_nd, c2t_st);
            // cout<<"s2c_dur: "<<s2c_dur<<" c2t_dur: "<<c2t_dur<<endl;
            if(curMinDist > s2c+c2t){
                curMinDist = s2c+c2t;
                newStartPos = j;
            }
            aps.pop_back();
        }
        // cout<<"curMinDist: "<<curMinDist<<" minDist: "<<minDist<<endl;
        //collect res.
        if(minDist >= curMinDist){
            minDist = curMinDist;
            aps.push_back(cand[newStartPos]);
        }
        if(minDist < curMinDist && i>0||newStartPos+1>=budget){
            for(int st=i; st<budget; st++){
                res50Dist[st+1].push_back(minDist);
                res50Time[st+1].push_back(1.0*MYMAX_INT_HALF);
            }
            break;
        }
        //next unchosen cand
        newStartPos++; 
        curMinDist = INF;
        auto apx_grdy_end = timerNd();
        double dur_time = dur(apx_grdy_end, apx_grdy_start);
        res50Dist[i+1].push_back(minDist);
        res50Time[i+1].push_back(dur_time);
    }
    cout<<"minDist apx: "<<minDist<<endl
    <<"ap size: "<<aps.size()<<endl
    <<"aps: "<<endl;
    for(auto ele:aps){
        cout<<ele<<" ";
    }
    cout<<endl;
}
void algo_apx_greedy_pre(vector<int>& terminal, vector<int>& cand, int budget, int idx_flag){
    res50Dist.resize(budget+1);
    res50Time.resize(budget+1);
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long ub = 0;
    for(auto ele:src){
        if(idx_flag){
            ub+=distanceQuery(ele, sink);
        }else{
            ub+=dij(ele, sink);
        }
    }
    cout<<"ub: "<<ub<<endl;
    int pos = 0;
    int newStartPos = 0;
    long long minDist = LLINF;
    long long curMinDist = LLINF;
    vector<int> aps;
    //In the search tree, we greedily adopt the pre-order traverse.
    //We start from the 1st rank cand until exhaustively traversed all given cands.
    //Or we stop when the next cand does not decrease the tot weight of the paths.
    for(int i=0; i<budget; i++){
        auto apx_grdy_start = timerSt();
        aps.push_back(cand[i]);
        vector<vector<pair<int, int>>> newG;
        unordered_map<int, vector<int>> umpii;
        long long s2c = assign_src2cand(terminal, aps, newG, umpii, idx_flag);
        long long c2t = prim(newG, newG.size());
        //collect res.
        if(minDist > s2c+c2t){
            curMinDist = s2c+c2t;
        }
        //it means no further comp is needed, 
        //as we can not further reduce the tot cost by increasing # of assembly points
        else
        {
            for (int st = i; st < budget; st++)
            {
                res50Dist[st + 1].push_back(minDist);
                res50Time[st + 1].push_back(1.0 * MYMAX_INT_HALF);
            }
            //no need for further comp, break
            break;
        }
        minDist = min(minDist, curMinDist);
        //reset for next unchosen cand
        curMinDist = INF;
        auto apx_grdy_end = timerNd();
        double dur_time = dur(apx_grdy_end, apx_grdy_start);
        res50Dist[i+1].push_back(minDist);
        res50Time[i+1].push_back(dur_time);
    }
    cout<<"minDist apx: "<<minDist<<endl
    <<"ap size: "<<aps.size()<<endl
    <<"aps: "<<endl;
}
void algo_greedy_split_O_n(vector<int>& terminal, vector<int>& cand, int budget){
    // auto time_st = timerSt();
    // vector<int> src(terminal.begin(), terminal.end()-1);
    // int sink = terminal.back();
    // //first for the number of srcs assigned to this cand/sink;
    // //second vector to store who are these srcs
    // //    the first the ele of the vector is the cand, and the rest are the srcs assigned
    // priority_queue<pair<int, pair<int, vector<int>>>> pq;
    // long long upBd = 0;
    // for(auto ele:src){
    //     upBd+=distanceQuery(ele, sink);
    // }
    // pq.emplace(make_pair(src.size(), make_pair(sink, src)));

    // vector<int> newCand;
    // int pos=0;
    // unordered_map<int, vector<int>> umpiv;
    // umpiv.emplace(make_pair(sink, src));
    // while(!pq.empty()){
    //     vector<int> uvec = pq.top().second.second;
    //     // the ap is not valueable if it has less than 2 srcs assigned
    //     if(uvec.size() < 2+1) break;
    //     int u = pq.top().second.first;
    //     pq.pop();
    //     //we relax the sink by cands and newly adopted cands
    //     if(u == sink){
    //         vector<int> assigned2Cand;
    //         //find the min cadn that minimize the tot most
    //         //@uvec we reuse this to return left srcs assigned to the root
    //         //@assigned2Cand we return the cands that achieves most decrease
    //         //@pos we need a new pos to greedily neglect visited cands.
    //         int min_cand = locate_min_cand_decrease_most(u, uvec, cand, newCand, assigned2Cand, pos);
    //         //no valuable cand, continue.
    //         if(-1 == min_cand) continue;
    //         // validate if spanning tree+greedy assignment can achieve a lower overall cost
    //         // @umpiv updated
    //         // @upBd updated
    //         bool ifVal = validate_overall_cost(min_cand, terminal, umpiv, upBd);
    //         if(!ifVal) continue; 
    //         //push to pq and will be processsed later
    //         pq.emplace(make_pair(umpiv[sink].size(), make_pair(sink, umpiv[sink])));
    //         if(umpiv.find(min_cand) != umpiv.end()){
    //             pq.emplace(make_pair(umpiv[min_cand].size(), make_pair(min_cand, umpiv[min_cand])));
    //         }
    //     }
    //     //deal with newly added cands and expand from the new cand to check if it 
    //     //has promising neighbors in the expansion.
    //     else{
    //         int min_cand = find_promising_nei_by_bfs(u, uvec, cand, newCand, budget);
    //         if(-1 == min_cand) continue;
    //         bool ifVal = validate_overall_cost(min_cand, terminal, umpiv, upBd);
    //         if(!ifVal) continue; 
    //     }
    // }
    // auto time_nd = timerNd();
    // double apx_time = dur(time_nd, time_st);
    // cout<<"Report:"<<endl;
    // cout<<"time: "<<apx_time<<endl;
    // cout<<"minDist: "<<upBd<<endl;
    // // for(auto ele:terminal){
    // //     cout<<ele<<" ";
    // // }
    // // cout<<endl;
    // for(auto ele:umpiv){
    //     cout<<ele.first<<" "<<ele.second.size()<<endl;
    // }
}
void algo_pure_split_5_8(vector<int>& terminal, vector<int>& cand, int budget, string whichFS){
    auto algo_st = timerSt();
    int kth_ap=0;
    vector<pair<int, pair<vector<int>, long long>>> term_assndSrc;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long minDist = 0;
    //CAUTION: vector<int> will cost too much memory when vector is large!!!
    //this pq use first ele of 1st pair as the rank
    //use vector<int>* to save memory
    priority_queue<pair<int, pair<int, pair<vector<int>, long long>>>> pq;
    //init tot
    long long tot=0;
    for(auto ele:src){
        tot += distanceQuery(sink, ele);
    }
    cout<<"upBd: "<<tot<<endl;
    //
    for(int loop=1; loop<=budget; loop++){
        while(!pq.empty()){
            pq.pop();
        }
        minDist = 0;
        term_assndSrc.clear();
        pq.emplace(make_pair(src.size(), make_pair(sink, make_pair(src, tot))));
        vector<int> vis(GraphSize+2, 0);
        // recursion formula: tot = splt1 + splt2 + e(splt1, splt2)
        while (!pq.empty())
        {
            int u = pq.top().second.first;
            int cnt = pq.top().first;
            long long tot = pq.top().second.second.second;
            vector<int> src2split = pq.top().second.second.first;
            // cout<<"*pop ver: "<<u<<" assigned num: "<<cnt<<" cur tot: "<<tot<<endl;
            if (cnt <= 2 || term_assndSrc.size() >= loop ||pq.size() > loop)
            {
                // collect all and break
                while (!pq.empty())
                {
                    int loc_u = pq.top().second.first;
                    int loc_cnt = pq.top().first;
                    long long loc_tot = pq.top().second.second.second;
                    vector<int> loc_src2split = pq.top().second.second.first;
                    pq.pop();
                    if (loc_tot == 0 && loc_u != sink)
                        continue;
                    cout<<"-> cand: "<<loc_u<<" cnt: "<<loc_cnt<<" tot: "<<loc_tot<<endl;
                    term_assndSrc.push_back(make_pair(loc_u, make_pair(loc_src2split, loc_tot)));
                }
                break;
            }
            pq.pop();

            vector<int> src2root;
            vector<int> src2cand;
            //tot is updated
            int min_cand = -1;
            if("bfs" == whichFS){
                min_cand = find_cand_bfs(u, tot, src2split);
            }
            if("dfs" == whichFS){
                min_cand = find_cand_dfs(u, tot, src2split);
            }
            if (-1 == min_cand)
            {
                // We check every cand and sink that can possibly be split srcs,
                // if it can no longer find a split plan, then we stop split on this ver.
                // the root poped has been fully relaxed. Store it.
                term_assndSrc.push_back(make_pair(u, make_pair(src2split, tot)));
            }
            else
            {
                vector<int> split2cand;
                vector<int> split2root;
                long long cand_tot = 0;
                long long root_tot = 0;
                // split and calc
                allocate_cand_root(u, min_cand, src2split, split2cand, split2root, cand_tot, root_tot);
                pq.emplace(make_pair(split2cand.size(), make_pair(min_cand, make_pair(split2cand, cand_tot))));
                pq.emplace(make_pair(split2root.size(), make_pair(u, make_pair(split2root, root_tot))));
            }
        }
        auto algo_nd = timerNd();
        double algo_exe_time = dur(algo_nd, algo_st);
        vector<vector<pair<int, int>>> newG;
        genG(newG, term_assndSrc);
        long long c2t = prim(newG, newG.size());
        cout<<"term size "<<term_assndSrc.size()<<endl;
        for(auto ele:term_assndSrc){
            minDist+=ele.second.second;
        }
        minDist+=c2t;
        cout<<loop<<"-minDist: "<<minDist<<" ap #: "<<term_assndSrc.size()-1<<" time "<<algo_exe_time<<endl;
        res50Time[loop].push_back(algo_exe_time);
        res50Dist[loop].push_back(minDist);
    }
    
}
void algo_pure_split_1e1_4(vector<int>& terminal, vector<int>& cand, int budget, string whichFS){
    auto algo_st = timerSt();
    int timeoutFlag = 0;
    int kth_ap=0;
    vector<pair<int, pair<vector<int>, long long>>> term_assndSrc;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long minDist = 0;
    priority_queue<pair<int, pair<int, pair<vector<int>, long long>>>> pq;
    //init tot
    long long tot=0;
    for(auto ele:src){
        tot += distanceQuery(sink, ele);
    }
    cout<<"upBd: "<<tot<<endl;
    //
    for(int loop=budget; loop<=budget; loop++){
        while(!pq.empty()){
            pq.pop();
        }
        vector<int> vis(GraphSize+2, 0);
        minDist = 0;
        term_assndSrc.clear();
        pq.emplace(make_pair(src.size(), make_pair(sink, make_pair(src, tot))));
        // recursion formula: tot = splt1 + splt2 + e(splt1, splt2)
        while (!pq.empty())
        {
            // timeoutFlag = timeOutCheck(algo_st);
            // if(1 == timeoutFlag) break;
            int u = pq.top().second.first;
            int cnt = pq.top().first;
            long long tot = pq.top().second.second.second;
            vector<int> src2split = pq.top().second.second.first;
            // cout<<"*pop ver: "<<u<<" assigned num: "<<cnt<<" cur tot: "<<tot<<endl;
            //term cond
            if (cnt <= 2 || term_assndSrc.size() >= loop ||pq.size() > loop)
            {
                // collect all and break
                while (!pq.empty())
                {
                    int loc_u = pq.top().second.first;
                    int loc_cnt = pq.top().first;
                    long long loc_tot = pq.top().second.second.second;
                    vector<int> loc_src2split = pq.top().second.second.first;
                    pq.pop();
                    if (loc_tot == 0 && loc_u != sink)
                        continue;
                    // cout<<"-> cand: "<<loc_u<<" cnt: "<<loc_cnt<<" tot: "<<loc_tot<<endl;
                    term_assndSrc.push_back(make_pair(loc_u, make_pair(loc_src2split, loc_tot)));
                }
                break;
            }
            pq.pop();

            vector<int> src2root;
            vector<int> src2cand;
            //@tot is updated
            // cout<<"√before: "<<tot<<endl;
            auto dfs_st = timerSt();
            int min_cand = -1;
            if("bfs" == whichFS){
                min_cand = find_cand_bfs(u, tot, src2split);
                if(min_cand == u){
                    min_cand = -1;
                }
                // auto bfs_nd = timerNd();
                // double bfs_tot = dur(bfs_nd, dfs_st);
                // cout<<"bfs_time: "<<bfs_tot<<endl;
            }
            if("dfs" == whichFS){
                min_cand = find_cand_dfs(u, tot, src2split);
                // auto dfs_nd = timerNd();
                // double dfs_tot = dur(dfs_nd, dfs_st);
                // cout<<"dfs_time: "<<dfs_tot<<endl;
            }
            // int min_cand = find_cand_dfs(u, tot, src2split, vis);
            // int min_cand = find_cand_bfs(u, tot, src2split, vis);
            // cout<<"×after: "<<tot<<endl;
            if (-1 == min_cand)
            {
                // We check every cand and sink that can possibly be split srcs,
                // if it can no longer find a split plan, then we stop split on this ver.
                // the root poped has been fully relaxed. Store it.
                term_assndSrc.push_back(make_pair(u, make_pair(src2split, tot)));
                
            }
            else
            {
                vector<int> split2cand;
                vector<int> split2root;
                long long cand_tot = 0;
                long long root_tot = 0;
                // split and calc
                allocate_cand_root(u, min_cand, src2split, split2cand, split2root, cand_tot, root_tot);
                // cout<<"s2cand: "<<split2cand.size()<<" s2root: "<<split2root.size()<<endl;
                pq.emplace(make_pair(split2cand.size(), make_pair(min_cand, make_pair(split2cand, cand_tot))));
                pq.emplace(make_pair(split2root.size(), make_pair(u, make_pair(split2root, root_tot))));
            }
        }
        if(0 == timeoutFlag){
            auto algo_nd = timerNd();
            double algo_exe_time = dur(algo_nd, algo_st);
            vector<vector<pair<int, int>>> newG;
            genG(newG, term_assndSrc);
            long long c2t = prim(newG, newG.size());
            cout << "term size " << term_assndSrc.size() << endl;
            for (auto ele : term_assndSrc)
            {
                // cout<<ele.first<<" assigned #: "<<ele.second.first.size()<<endl;
                minDist += ele.second.second;
            }
            minDist += c2t;
            cout << loop << "-minDist: " << minDist << " ap #: " << term_assndSrc.size() - 1 << " time " << algo_exe_time << endl;

            res50Time[loop].push_back(algo_exe_time);
            res50Dist[loop].push_back(minDist);
        }else{
            //timeout
            res50Time[loop].push_back(7200.00);
            res50Dist[loop].push_back(LLONG_MAX);
        }
    }
    
}
void run_AdvDP(int code, char* file, int k){
    //read query file
    FILE* pQf = fopen(file, "r");
    vector<int> Q;
    int tmpNodeID = 0;
    vector<vector<int>> terminals;
    loadTerminals(terminals, file, code);
        
    for(int i=0; i<terminals.size(); i++){
        cout<<"---------------"<<i<<"--------------"<<endl;;
        algo6_AdvDP(terminals[i], k);
        cout<<"____----_____----____----____----"<<endl;
    }
    printRes2File();
}
void run_AdvDP_given_candidate(int code, char* file, char* cand_file, int k){
    //read query file
    FILE* pQf = fopen(file, "r");
    vector<int> Q;
    int tmpNodeID = 0;
    vector<vector<int>> terminals;
    while(EOF != fscanf(pQf, "%d", &tmpNodeID)){
        Q.push_back(tmpNodeID);
    }
    fclose(pQf);
    vector<int>  tmp;
    for(int j=0; j<(int)Q.size(); j++){
        tmp.push_back(Q[j]);
        if((int)tmp.size() == code){
            terminals.push_back(tmp);
            tmp.clear();
            continue;
        }
    }
    // for(auto line:terminals){
    //     for(auto ele:line){
    //         cout<<ele<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    //read candidate note that multiple whitespace between two cand vertices will cause error 
    fstream fs_cand(cand_file, ios::in);
    vector<unordered_set<int>> candidates; 
    string line;
    while(getline(fs_cand, line)){
        // cout<<"candline "<<line<<endl;
        stringstream ss(line);
        string cand;
        unordered_set<int> usi;
        while(getline(ss, cand, ' ')){
            usi.emplace(stoi(cand));
        }
        candidates.push_back(usi);
    }
    fs_cand.close();
    // cout<<"--cand--"<<endl;
    // for(auto ele:candidates[0]){
    //     cout<<ele<<" ";
    // }
    // cout<<endl;
    for(int i=0; i<terminals.size(); i++){
        cout<<"---------------"<<i<<"---------------"<<endl;;
        algo6_AdvDP_given_candidates(candidates[i], terminals[i], k);
        cout<<"____----____----____----____----"<<endl;
    }
    printRes2File();
}
void run_BasicDP(int code, char* file, int k){
    //read query file
    vector<vector<int>> terminals;
    loadTerminals(terminals, file, code);
        
    for(int i=0; i<terminals.size(); i++){
        cout<<"---------------"<<i<<"--------------"<<endl;;
        algo5_BasicDP(terminals[i], k);
        cout<<"____----_____----____----____----"<<endl;
    }
    printRes2File();
}


void run_BnB(int code, char* file, int k){
    //read query file
    vector<vector<int>> terminals;
    loadTerminals(terminals, file, code); 

    for(int i=0; i<terminals.size(); i++){
        cout<<"---------------"<<i<<"--------------"<<endl;;
        algo4_BnB(terminals[i], k, 0);
        cout<<"____----_____----____----____----"<<endl;
    }
    printRes2File();
}
void run_BnB_with_index(int code, char* file, int k){
    //read query file
    vector<vector<int>> terminals;
    loadTerminals(terminals, file, code); 

    for(int i=0; i<terminals.size(); i++){
        cout<<"---------------"<<i<<"--------------"<<endl;;
        algo4_Branch_and_Bound(terminals[i], k);
        algo4_BnB(terminals[i], k ,1);
        cout<<"____----_____----____----____----"<<endl;
    }
    printRes2File();
}

//bfs:sim dfs:simapp
void run_appx(int srcNsink_num, char* queryFile, string icf, int budget){
    res50Dist.resize(budget+1);
    res50Time.resize(budget+1);
    //read query file to terminals
    vector<vector<int>> terminals;
    
    loadTerminals(terminals, queryFile, srcNsink_num);
    //gen cands
    vector<vector<int>> cand;
    cand.resize(terminals.size()+1);
    vector<string> apprxAlgo={"bfs","dfs"};
    string whichFS;
    for(int i=0; i<apprxAlgo.size();i++){
        if(icf.find(apprxAlgo[i]) != string::npos){
            whichFS = apprxAlgo[i];
        }
    }
    for(int i=0; i<terminals.size(); i++){
        auto tot_start = timerSt(); 

        //algo1 is not compatible to appx.
        auto algo1_st = timerSt(); 
        cout<<"---------------"<<i<<"--------------"<<endl;;
        // algo_apx_greedy(terminals[i], cand[i], budget);
        // algo_greedy_split_O_n(terminals[i], cand[i], budget);
        // algo_apx_greedy_pre(terminals[i], cand[i], budget);
        if(budget < 10){
            algo_pure_split_5_8(terminals[i], cand[i], budget, whichFS);
        }else{
            algo_pure_split_1e1_4(terminals[i], cand[i], budget, whichFS);
        }
        cout<<"____----_____----____----____----"<<endl;
        auto algo1_nd = timerNd(); 
        double tot_algo1 = dur(algo1_nd, algo1_st);
        
        auto tot_end = timerNd(); 
        double tot_time = dur(tot_end, tot_start);
        cout\
        <<" tot_idx: "<<tot_idx_time\
        <<" tot_path_time: "<<tot_get_paths_time\
        <<" tot_other_1 : "<<tot_other_1\
        <<" tot_other_2 : "<<tot_other_2\
        <<endl

        <<" tot_algo1: "<<tot_algo1\
        <<endl

        <<" tot_time: "<<tot_time\
        <<endl;
        cout<<endl;
        tot_idx_time = 0.0;
        tot_get_paths_time = 0.0;
    }
    printRes2File();
}

void run_apx_cand(int srcNsink_num, char* queryFile, string icf, int budget, int idx_flag){
    res50Dist.resize(budget+1);
    res50Time.resize(budget+1);
    //read query file to terminals
    vector<vector<int>> terminals;
    loadTerminals(terminals, queryFile, srcNsink_num);
    //gen cands
    vector<vector<int>> cand;
    cand.resize(terminals.size()+1);
    string bet = "bet";
    string pcp = "pcp";
    string mst = "mst"; 
    for(int i=0; i<terminals.size(); i++){
        //gen cand
        auto tot_start = timerSt(); 

        auto select_cand_start = timerSt(); 
        //bet sorted by # of sps passes through each vertex
        if(icf.find(bet) != string::npos){
            if(idx_flag){
                gen_bet_with_idx(terminals[i], cand[i]);
            }else{
                gen_bet(terminals[i], cand[i]);
            }
        }
        //ver on pcp ranked by degree
        if(icf.find(pcp) != string::npos){
            if(idx_flag){
                gen_pcp_with_idx(terminals[i], cand[i]);
            }else{
                gen_pcp(terminals[i], cand[i]);
            }
        }
        //ver on mst ranked by degree
        if(icf.find(mst) != string::npos){
            gen_mst(terminals[i], cand[i]);
        }
        auto select_cand_end = timerNd(); 
        double tot_select_cand = dur(select_cand_end, select_cand_start); 

        //algo1 is not compatible to appx.
        //now we have cands
        //we simply find its cost 
        auto a_st= timerSt(); 
        cout<<"---------------"<<i<<"--------------"<<endl;;
        algo_apx_greedy_pre(terminals[i], cand[i], budget, idx_flag);
        cout<<"____----_____----____----____----"<<endl;
        auto a_nd = timerNd(); 
        double tot_a = dur(a_nd, a_st);
        
        auto tot_end = timerNd(); 
        double tot_time = dur(tot_end, tot_start);
        cout\
        <<" select cand time: "<<tot_select_cand\

        <<" tot_algo: "<<tot_a\
        <<endl

        <<" tot_time: "<<tot_time\
        <<endl;
        cout<<endl;
        //refresh
        tot_idx_time = 0.0;
        tot_get_paths_time = 0.0;
    }
    //store res to file
    printRes2File();
    //store cands to file
    printCand2File(icf, cand);
}
void run_simpleGreedy(int srcNsink_num, char* queryFile, string icf, int budget){
    // vector<int>& terminal, int budget, vector<int>& aps, queue<int>& qi, int& upDist
    vector<vector<int>> terminals;
    loadTerminals(terminals, queryFile, srcNsink_num);
    vector<vector<int>> APset; 
    APset.resize(terminals.size()+1);

    for(int i=0; i<terminals.size(); i++){
        vector<int> srcs(terminals[i].begin(), terminals[i].end()-1);
        long long upbd = terminals_cand_tot_dis(srcs, terminals[i].back());
        auto res = algo_simpleGreedy(terminals[i], budget, APset[i], upbd);
        cout<<"res of terminals["<<i<<"]: "<<res<<endl;
    }
}
// How to run
//            |----binary----|-----------------------------------------input------------------------------------------------------|-------------output-------------|
// meaning    |exe           |index       |graph      |#(src+sink)     |query        |budget     |algo           |candidates(i/o) |res 2 file    |redirect to log  |    
// var in code|              |infileIDX   |infileGR   |srcNsink_num    |queryFile    |budget     |chooseAlgo     |candFile        |outfile       |                 |
// input args |              |argv[1]     |argv[2]    |argv[3]         |argv[4]      |argv[5]    |argv[6]        |argv[7]         |argv[8]       |                 |
// shell       ./MSPA.out     NY.index     NY.gr        6               NY6.Q         5           3               NY6_[].cand      NY6.res        > NY6.log        |
int main(int argc, char** argv){

    //read index and graph
    char* infileIDX = argv[1];
    char* infileGR = argv[2];
    string ifidx = infileIDX;
    string ifgr = infileGR;
    int chooseAlgo = 0;
    chooseAlgo = char2int(argv[6]);
    if(1 == chooseAlgo || 2==chooseAlgo || 5==chooseAlgo || 6==chooseAlgo){
        readIndex(infileIDX);
    }
	readGraph(infileGR);
    cout<<"index: "<<ifidx<<endl;
    cout<<"graph: "<<ifgr<<endl;

    //params for MSPA
    int srcNsink_num = char2int(argv[3]);
    cout<<"#(src+sink): "<<srcNsink_num<<endl;

    char* queryFile = argv[4];
    string qf = queryFile;
    cout<<"query: "<<qf<<endl;

    int budget = 0; 
    budget = char2int(argv[5]);
    cout<<"budget: "<<budget <<endl;
    
    vector<string> algosVec={"", "BnB", "BasicDP", "AdvDP", "AdvDP given cands", \
    "bfs, dfs", "bet, pcp, mst", "naive greedy"};
    cout<<"algo: "<<algosVec[chooseAlgo]<<endl;
    
    //works only for "4 adv dp given cands"
    incandfile = argv[7];
    string icf = incandfile;

    outfile = argv[8];
    string of = outfile;
    outOpenFile = fopen(outfile, "w");
    
    cout<<"candidates: "<<icf<<endl;
    cout<<"res file: "<<of<<endl;
    
    // run algos

    // 01: branch and bound, index free or with index
    // we use int in this algo, as we only conduct on graphs with much less vertices
    if(0 == chooseAlgo){
        run_BnB(srcNsink_num, queryFile, budget);
    }
    if(1 == chooseAlgo){
        run_BnB_with_index(srcNsink_num, queryFile, budget);
    }
    // 2: Basic DP requires index.
    if(2 == chooseAlgo){
        run_BasicDP(srcNsink_num, queryFile, budget);
    }
    // 3: Adv DP is index free
    if(3 == chooseAlgo){
        run_AdvDP(srcNsink_num, queryFile, budget);
    }
    // 4:Run Adv DP given candidates
    if(4 == chooseAlgo){ 
        run_AdvDP_given_candidate(srcNsink_num, queryFile, incandfile, budget);
    }
    // 5: bfs and dfs based apprx algos. For 
    if(5 == chooseAlgo){
        run_appx(srcNsink_num, queryFile, icf, budget);
    }
    // 6: gen cands for mst, pcp and bet
    // we stop computation if more cands causes more than less cands.
    // as the cands for mst, pcp and bet are already ranked, we pick ap set sequentially
    // e.g {v1, v2, v3, v4, v5}, we pick {v1}, {v1, v2}, {v1, v2, v3}, ... as budget for 1, 2, 3, ...
    if(6 == chooseAlgo){
        run_apx_cand(srcNsink_num, queryFile, icf, budget, 1);
    }
    // 7: it is based on a naive idea that each newly selected assembly point should 
    // relax one of the paths connected  to the already selected assembly point
    // TODO:
    if(7 == chooseAlgo){
        run_simpleGreedy(srcNsink_num, queryFile, icf, budget);
    }
    fclose(outOpenFile);
    return 0;
}