#include "MSPAUtilities.h"
#include "pruned_landmark_labeling.h"
#include "H2H.h"
#include "appx.h"
using namespace std;

int algo1_comp_all_path(vector<int>& apset, vector<int>& terminal, vector<pair<int, long long>>& Gprim){
    int tot_dis=0;
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    int s2c = 0;
    int c2t = 0;
    for(int i=0; i<srcs.size(); i++){
        int min_Src2cand = INT_MAX; 
        for(int j=0; j<apset.size(); j++){
            int src2cand = distanceQuery(apset[j], srcs[i]);
            if(min_Src2cand > src2cand){
                min_Src2cand = src2cand;
            }
        }
        s2c+=min_Src2cand;
    }
    apset.push_back(sink);
    vector<vector<pair<int, int>>> newG;
    fillG(newG, apset);
    c2t = prim(newG, apset.size());
    apset.pop_back();
    tot_dis += s2c;
    tot_dis += c2t;
    return tot_dis;
}

void search(vector<int>& apset, int k, vector<pair<int, long long>>& Gprim, vector<int>& terminal, int pos, int& upBd, vector<int>& rstar){
    int given_ap_tot_dis = algo1_comp_all_path(apset, terminal, Gprim);
    if(apset.size()<k){
        for(int i=pos; i<GraphSize; i++){
            apset.push_back(Gprim[i].first);
            if(lowerBound(apset, terminal) < given_ap_tot_dis){
                search(apset, k, Gprim, terminal, i+1, upBd, rstar);
            }
            apset.pop_back();
        }
    }else if(apset.size() == k){
        if(given_ap_tot_dis < upBd){
            upBd = given_ap_tot_dis;
            cout<<upBd<<endl;
            rstar = apset;
            for(auto ele:rstar){
                cout<<ele<<" ";
            }
            cout<<endl;
        }
    }
}

void algo4_Branch_and_Bound(vector<int>& terminals, int k){
    //reorder by dis
    vector<pair<int, long long>> reorder_by_dis;
    for(int i=1; i<=GraphSize; i++){
        reorder_by_dis.push_back({i, terminals_cand_tot_dis(terminals, i)});
    }
    sort(reorder_by_dis.begin(), reorder_by_dis.end(), [&](pair<int, int> pa, pair<int, int> pb){
        return pa.second<pb.second;
    });
    vector<int> apset;
    int upBd = 0;
    int sink = terminals.back();
    vector<int> srcs(terminals.begin(), terminals.end()-1);
    for(auto ele:srcs){
        upBd += distanceQuery(ele, sink);
    }
    vector<int> rstar;
    // apset.push_back(sink);
    search(apset, k, reorder_by_dis, terminals, 0, upBd, rstar);
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
                    dpQ[i][j][kth][l] = TRIQUATINF;
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
    }
    
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
                dp[i][j][kth] = TRIQUATINF;
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

void algo_apx_greedy_pre(vector<int>& terminal, vector<int>& cand, int budget){
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
    //In the search tree, we greedily adopt the pre-order traverse.
    //We start from the 1st rank cand until exhaustively traversed all given cands.
    //Or we stop when the next cand does not decrease the tot weight of the paths.
    for(int i=0; i<budget; i++){
        auto apx_grdy_start = timerSt();
        aps.push_back(cand[i]);
        vector<vector<pair<int, int>>> newG;
        // auto s2c_st = timerSt();
        unordered_map<int, vector<int>> umpii;
        long long s2c = assign_src2cand(terminal, aps, newG, umpii);
        // auto s2c_nd = timerNd();
        // auto c2t_st = timerSt();
        long long c2t = prim(newG, newG.size());
        // auto c2t_nd = timerNd();
        // double s2c_dur = dur(s2c_nd, s2c_st);
        // double c2t_dur = dur(c2t_nd, c2t_st);
        // cout<<"s2c_dur: "<<s2c_dur<<" c2t_dur: "<<c2t_dur<<endl;
        // cout<<"curMinDist: "<<curMinDist<<" minDist: "<<minDist<<endl;
        //collect res.
        if(minDist > s2c+c2t){
            curMinDist = s2c+c2t;
        }else{
            auto pre_end = timerNd();
            double pre_dur_time = dur(pre_end, apx_grdy_start);
            res50Dist[i+1].push_back(minDist);
            res50Time[i+1].push_back(pre_dur_time);
            break;
        }
        minDist = min(minDist, curMinDist);
        //next unchosen cand
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
// we grdily split heavily load root with bfs/dfs strategy
void algo_pure_split_5_8(vector<int>& terminal, vector<int>& cand, int budget, string whichFS){
    int kth_ap=0;
    // collection of assembly points
    // p1 1st: root
    // p2 1st: assigned src set
    // @update 23.12.04 p3 1st: prev root
    // p2 2nd: subtree weight
    // vector<pair<int, pair<vector<int>, pair<int, long long>>>> term_assndSrc;
    vector<Spliter> term_assndSrc;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long minDist = 0;
    // p1 1st: num of src assigned to curRoot
    // p2 1st: curRoot
    // p3 1st: src set assgned to curRoot
    // @update 23.12.04 p4 1st: remember its prev root
    // p4 2nd: tot subtree weight rooted at curRoot
    // priority_queue<pair<int, pair<int, pair<vector<int>, pair<int, long long>>>>> pq;
    priority_queue<pair<int, Spliter>, vector<pair<int, Spliter>>, ComparePairs<Spliter>> pq;
    //init 1st ele enque
    long long tot=0;
    for(auto ele:src){
        tot += distanceQuery(sink, ele);
    }
    cout<<"upBd: "<<tot<<endl;
    // Given budget, we compute all possible cands within the budget limit.
    // Therefore, we use a loop for budget=1~5
    for(int loop=1; loop<=budget; loop++){
        long long loopTot = tot;
        //init
        auto algo_st = timerSt();
        while(!pq.empty()){
            pq.pop();
        }
        minDist = 0;
        term_assndSrc.clear();
        Spliter rs;
        rs.assignSrcSet = src;
        rs.subtreeRoot = sink;
        rs.prevRoot = sink;
        rs.subtreeWeight = loopTot;
        pq.emplace(make_pair(src.size(), rs));
        vector<int> vis(GraphSize+2, 0);
        // recursion formula: tot = splt1 + splt2 + e(splt1, splt2)
        while (!pq.empty())
        {
            // int u = pq.top().second.first;
            // int cnt = pq.top().first;
            // long long tot = pq.top().second.second.second.second;
            // int prev = pq.top().second.second.second.first;
            // vector<int> src2split = pq.top().second.second.first;
            Spliter spliter = pq.top().second;
            int u = spliter.subtreeRoot;
            int cnt = pq.top().first;
            long long subTot = spliter.subtreeWeight;
            int prev = spliter.prevRoot;
            vector<int> src2split = spliter.assignSrcSet;
            // limit reached and break
            // |S| <= k; budget limit;
            if (cnt <= 2 || term_assndSrc.size() >= loop ||pq.size() > loop)
            {
                // collect all and break
                while (!pq.empty())
                {
                    // int loc_u = pq.top().second.first;
                    // int loc_cnt = pq.top().first;
                    // long long loc_tot = pq.top().second.second.second;
                    // vector<int> loc_src2split = pq.top().second.second.first;
                    Spliter loc_rs = pq.top().second;
                    int loc_u = loc_rs.subtreeRoot;
                    int loc_cnt = pq.top().first;
                    long long loc_tot = loc_rs.subtreeWeight;
                    vector<int> loc_src2split = loc_rs.assignSrcSet;
                    pq.pop();
                    if (loc_tot == 0 && loc_u != sink)
                        continue;
                    cout<<"-> cand: "<<loc_u<<" cnt: "<<loc_cnt<<" tot: "<<loc_tot<<endl;
                    // term_assndSrc.push_back(make_pair(loc_u, make_pair(loc_src2split, loc_tot)));
                    term_assndSrc.push_back(loc_rs);
                }
                break;
            }
            pq.pop();

            vector<int> src2root;
            vector<int> src2cand;
            //@tot is updated
            //In App and SimApp we split srcs to root and best cand found.
            //The tree rooted at the root is then divided into 3 parts, subtree 1: src set 1 to cand, subtree2: src set 2 to root, and dist(cand, root).
            //In the next round, we grdily use the weight of the weight of a tree rooted at cand(or root) to relax the weight of subtree1 (or subt2).
            int min_cand = -1;
            if("bfs" == whichFS){
                min_cand = find_cand_bfs(u, loopTot, src2split);
                //TODO: test if init_fs and gen exclude will help accelerate and avoid scanning whole graph is necessary
                if(src.size() == src2split.size() && -1 == min_cand){
                    // vector<int> exclude(GraphSize+1, 0);
                    // min_cand = init_fs(u, loopTot, src2split, exclude); 
                    min_cand = init_fs(u, loopTot, src2split); 
                }
                    // min_cand = find_cand_bfs(u, loopTot, src2split);
            }
            // only vis once is bad strategy especially for bfs
            // int min_cand = find_cand_bfs(u, tot, src2split, vis);
            if("dfs" == whichFS){
                min_cand = find_cand_dfs(u, loopTot, src2split);
                //TODO: test if init_fs and gen exclude will help accelerate and avoid scanning whole graph is necesary
                if(src.size() == src2split.size() && -1 == min_cand){
                    // vector<int> exclude(GraphSize+1, 0);
                    // min_cand = init_fs(u, loopTot, src2split, exclude); 
                    min_cand = init_fs(u, loopTot, src2split); 
                }
            }
            // int min_cand = find_cand_dfs(u, tot, src2split, vis);
            if (-1 == min_cand)
            {
                // We check every cand and sink that can possibly be split srcs,
                // if it can no longer find a split plan, then we stop split on this ver.
                // By now, the root poped has been fully relaxed. Store it.
                // term_assndSrc.push_back(make_pair(u, make_pair(src2split, tot)));
                term_assndSrc.push_back(spliter);
                cout<<"-> cand: "<<u<<" cnt: "<<cnt<<" tot: "<<subTot<<endl;
            }
            else
            {
                vector<int> split2cand;
                vector<int> split2root;
                long long cand_tot = 0;
                long long root_tot = 0;
                // split and calc
                // It is possible that all src assigned to a new root,
                // the heavy loaded new root will be relaxed in the future.
                allocate_cand_root(u, min_cand, src2split, split2cand, split2root, cand_tot, root_tot);
                Spliter cand_rs(split2cand, min_cand, u, cand_tot);
                Spliter root_rs(split2root, u, prev, root_tot);
                pq.emplace(make_pair(split2cand.size(), cand_rs));
                pq.emplace(make_pair(split2root.size(), root_rs));
            }
        }
        // @update 23.12.04 due to prev root record, prim is no longer needed.
        // vector<vector<pair<int, int>>> newG;
        // genG(newG, term_assndSrc);
        // long long c2t = prim(newG, newG.size());
        // cout<<"term size "<<term_assndSrc.size()<<endl;
        for(auto ele:term_assndSrc){
            // cout<<ele.first<<" assigned #: "<<ele.second.first.size()<<endl;
            minDist+=ele.subtreeWeight;
            minDist+=distanceQuery(ele.subtreeRoot, ele.prevRoot);
        }
        auto algo_nd = timerNd();
        // minDist+=c2t;
        double algo_exe_time = dur(algo_nd, algo_st);
        cout<<loop<<"-minDist: "<<minDist<<" ap #: "<<term_assndSrc.size()-1<<" time "<<algo_exe_time<<endl;
        res50Time[loop].push_back(algo_exe_time);
        res50Dist[loop].push_back(minDist);
    }
    
}
void algo_pure_split_1e1_4(vector<int>& terminal, vector<int>& cand, int budget, string whichFS){
    auto algo_st = timerSt();
    int timeoutFlag = 0;
    int kth_ap=0;
    // vector<pair<int, pair<vector<int>, long long>>> term_assndSrc;
    vector<Spliter> term_assndSrc;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    long long minDist = 0;
    //@updated 23.12.04 struct to replace nested pairs.
    // priority_queue<pair<int, pair<int, pair<vector<int>, long long>>>> pq;
    priority_queue<pair<int, Spliter>, vector<pair<int, Spliter>>, ComparePairs<Spliter>> pq;
    //init tot
    long long tot=0;
    for(auto ele:src){
        tot += distanceQuery(sink, ele);
    }
    cout<<"upBd: "<<tot<<endl;
    //the budget is fixed, so no more loops
    for(int loop=budget; loop<=budget; loop++){
        while(!pq.empty()){
            pq.pop();
        }
        vector<int> vis(GraphSize+2, 0);
        minDist = 0;
        term_assndSrc.clear();
        Spliter rs(src, sink, sink, tot);
        // pq.emplace(make_pair(src.size(), make_pair(sink, make_pair(src, tot))));
        pq.emplace(make_pair(src.size(), rs));
        // recursion formula: tot = splt1 + splt2 + e(splt1, splt2)
        while (!pq.empty())
        {
            timeoutFlag = timeOutCheck(algo_st, 1800);
            // if(1 == timeoutFlag) break;
            // int u = pq.top().second.first;
            // int cnt = pq.top().first;
            // long long tot = pq.top().second.second.second;
            // vector<int> src2split = pq.top().second.second.first;
            // cout<<"*pop ver: "<<u<<" assigned num: "<<cnt<<" cur tot: "<<tot<<endl;
            int cnt = pq.top().first;
            Spliter rs = pq.top().second;
            int u = rs.subtreeRoot;
            long long tot = rs.subtreeWeight;
            vector<int> src2split = rs.assignSrcSet;
            int prev = rs.prevRoot;

            //termination condition
            //If timeout, we keep the cur found aps, and return;
            if (timeoutFlag||cnt <= 2 || term_assndSrc.size() >= loop ||pq.size() > loop)
            {
                // collect all and break
                while (!pq.empty())
                {
                    // int loc_u = pq.top().second.first;
                    // int loc_cnt = pq.top().first;
                    // long long loc_tot = pq.top().second.second.second;
                    // vector<int> loc_src2split = pq.top().second.second.first;
                    Spliter  loc_rs = pq.top().second;
                    long long loc_tot = loc_rs.subtreeWeight;
                    int loc_u = loc_rs.subtreeRoot;
                    pq.pop();
                    if (loc_tot == 0 && loc_u != sink)
                        continue;
                    // cout<<"-> cand: "<<loc_u<<" cnt: "<<loc_cnt<<" tot: "<<loc_tot<<endl;
                    // term_assndSrc.push_back(make_pair(loc_u, make_pair(loc_src2split, loc_tot)));
                    term_assndSrc.push_back(loc_rs);
                }
                break;
            }
            pq.pop();

            vector<int> src2root;
            vector<int> src2cand;
            //@tot is updated
            auto dfs_st = timerSt();
            int min_cand = -1;
            if("bfs" == whichFS){
                min_cand = find_cand_bfs(u, tot, src2split);
            }
            if("dfs" == whichFS){
                min_cand = find_cand_dfs(u, tot, src2split);
            }
            // Use vis to vis each ver only once. Much faster, but worse.
            // int min_cand = find_cand_dfs(u, tot, src2split, vis);
            // int min_cand = find_cand_bfs(u, tot, src2split, vis);
            if (-1 == min_cand)
            {
                // We check every cand and sink that can possibly be split srcs,
                // if it can no longer find a split plan, then we stop split on this ver.
                // the root poped has been fully relaxed. Store it.
                // term_assndSrc.push_back(make_pair(u, make_pair(src2split, tot)));
                   term_assndSrc.push_back(rs); 
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
                // pq.emplace(make_pair(split2cand.size(), make_pair(min_cand, make_pair(split2cand, cand_tot))));
                // pq.emplace(make_pair(split2root.size(), make_pair(u, make_pair(split2root, root_tot))));
                Spliter cand_rs(split2cand, min_cand, u, cand_tot);
                Spliter root_rs(split2root, u, prev, root_tot);
                pq.emplace(make_pair(split2cand.size(), cand_rs));
                pq.emplace(make_pair(split2root.size(), root_rs));
            }
        }
        // if(0 == timeoutFlag){
            auto algo_nd = timerNd();
            double algo_exe_time = dur(algo_nd, algo_st);
            // @updated 23.12.04 prev root recorded, prim no longer needed.
            // vector<vector<pair<int, int>>> newG;
            // genG(newG, term_assndSrc);
            // long long c2t = prim(newG, newG.size());
            // cout << "term size " << term_assndSrc.size() << endl;
            // minDist += c2t;
            for (auto ele : term_assndSrc)
            {
                // cout<<ele.first<<" assigned #: "<<ele.second.first.size()<<endl;
                minDist += ele.subtreeWeight;
                minDist += distanceQuery(ele.prevRoot, ele.subtreeRoot);
            }
            //term_assndSrc contains the sink
            cout << loop << "-minDist: " << minDist << " ap #: " << term_assndSrc.size() - 1 << " time " << algo_exe_time << endl;

            res50Time[loop].push_back(algo_exe_time);
            res50Dist[loop].push_back(minDist);
        // }
        // else{
        //     // timeout
        //     res50Time[loop].push_back(7200.00);
        //     res50Dist[loop].push_back(LLONG_MAX);
        // }
    }
    
}
void run_AdvDP(int code, char* file, int k){
    //read query file
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
    vector<vector<int>> terminals;
    loadTerminals(terminals, file, code);
    
    //read candidate note that multiple whitespace between two cand vertices will cause error 
    fstream fs_cand(cand_file, ios::in);
    vector<unordered_set<int>> candidates; 
    string line;
    while(getline(fs_cand, line)){
        stringstream ss(line);
        string cand;
        unordered_set<int> usi;
        while(getline(ss, cand, ' ')){
            usi.emplace(stoi(cand));
        }
        candidates.push_back(usi);
    }
    fs_cand.close();
    
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
        algo4_Branch_and_Bound(terminals[i], k);
        cout<<"____----_____----____----____----"<<endl;
    }
    printRes2File();
}

void run_appx(int srcNsink_num, char* queryFile, string icf, int budget){
    res50Dist.resize(budget+1);
    res50Time.resize(budget+1);
    //read query file to terminals
    vector<vector<int>> terminals;
    loadTerminals(terminals, queryFile, srcNsink_num);
    vector<vector<int>> cand;
    cand.resize(terminals.size()+1);
    string bfs = "bfs";
    string dfs = "dfs";
    string whichFS;
    if(icf.find(bfs) != string::npos){
        whichFS=bfs;
    }
    if(icf.find(dfs) != string::npos){
        whichFS=dfs;
    }
    for(int i=0; i<terminals.size(); i++){
        auto tot_start = timerSt(); 
        cout<<"---------------"<<i<<"--------------"<<endl;;
        if(budget < 10){
            algo_pure_split_5_8(terminals[i], cand[i], budget, whichFS);
        }else{
            algo_pure_split_1e1_4(terminals[i], cand[i], budget, whichFS);
        }
        cout<<"____----_____----____----____----"<<endl;
        auto tot_end = timerNd(); 
        double tot_time = dur(tot_end, tot_start);
        cout<<i<<"th "<<whichFS<<" split time: "<<tot_time<<endl;
    }
    //store res to file
    printRes2File();
    //store cands to file
    printCand2File(icf, cand);
}

void run_apx_cand(int srcNsink_num, char* queryFile, string icf, int budget){
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
        cout<<"---------------"<<i<<"--------------"<<endl;;
        cout<<"top 10 ranked [id, val]: "<<endl;
        if(icf.find(bet) != string::npos){
            gen_bet_with_idx(terminals[i], cand[i]);
        }
        if(icf.find(pcp) != string::npos){
            gen_pcp_with_idx(terminals[i], cand[i]);
        }
        if(icf.find(mst) != string::npos){
            gen_mst(terminals[i], cand[i]);
        }
        auto select_cand_end = timerNd(); 
        double tot_select_cand = dur(select_cand_end, select_cand_start); 

        auto algo_apx_cand_st = timerSt(); 
        // Since the cands are given and ranked, we grdily traverse along the search tree in a pre-traversal order.
        // We traverse along the left most search tree.
        // If the tot weigth no longer decreases, we just break and retrun.
        algo_apx_greedy_pre(terminals[i], cand[i], budget);
        cout<<"____----_____----____----____----"<<endl;
        auto algo_apx_cand_nd = timerNd(); 
        double tot_algo1 = dur(algo_apx_cand_nd , algo_apx_cand_st );
        
        auto tot_end = timerNd(); 
        double tot_time = dur(tot_end, tot_start);
        cout<<"slct cand time: "<<tot_select_cand \
            <<", apx comp time: "<<tot_algo1\
            <<", tot time: "<<tot_time\
        <<endl;
    }
    //store res to file
    printRes2File();
    //store cands to file
    printCand2File(icf, cand);
}
// @param chooseAlgo
//        1ca: BnB; 2ca: BasicDP; 3ca: AdvDP; 4ca: AdvDP given Candidates; 
//        5ca: App(BFS-based) and SimApp(DFS-based) 6ca: grdy while comp cands online
// @param incandfile
//        useless for algo 1ca~3ca
//        need to be precomputed for algo 4ca
//        output cands computed during exe of algo 6ca
// @param infileIDX
//        need to be precomputed by tools provdide in folder tools
//        for pll, bit parallel num should be set as 64
// @param infileGR
//        need to be in format:
//     ```
//        n m                 //num of nodes and edges at the first line
//        u_1 v_1 w_1         //edge from u_1 to v_1, weight w_1. If social network, it should be 1
//        u_2 v_2 w_2         //for all graphs, nodes id starts from 1
//        ...
//        u_n v_n w_n
//     ``` 
// @param queryFile
//        precomputed with tool provided in folder tools
// @param outfile
//        results are stored here
// @param log
//        in the bash script, stdout and stderr are redirected to the log files.
//
// How to run
//            |----binary----|-----------------------------------------input------------------------------------------------------|-------------output-------------|
// meaning    |exe           |index       |graph      |#(src+sink)     |query        |budget     |algo           |candidates(i/o) |res 2 file    |redirect to log  |    
// var in code|              |infileIDX   |infileGR   |srcNsink_num    |queryFile    |budget     |chooseAlgo     |candFile        |outfile       |                 |
// input args |              |argv[1]     |argv[2]    |argv[3]         |argv[4]      |argv[5]    |argv[6]        |argv[7]         |argv[8]       |                 |
// shell       ./MSPA.out     NY.index     NY.gr        6               NY6.Q         5           3               NY6_[].cand      NY6.res        > NY6.log        |
int main(int argc, char** argv){
    char* infileIDX = argv[1];
    char* infileGR = argv[2];
    string ifidx = infileIDX;
    string ifgr = infileGR;
    int srcNsink_num = char2int(argv[3]);
    char* queryFile = argv[4];
    string qf = queryFile;
    int budget = 0; 
    budget = char2int(argv[5]);
    int chooseAlgo = 0;
    chooseAlgo = char2int(argv[6]);
    incandfile = argv[7];
    outfile = argv[8];
    string icf = incandfile;
    string of = outfile;

    cout<<"index file: "<<ifidx<<endl;
    cout<<"graph file: "<<ifgr<<endl;
    cout<<"#(src+sink): "<<srcNsink_num<<endl;
    cout<<"query file: "<<qf<<endl;
    cout<<"#budget: "<<budget <<endl;
    cout<<"algo: "<<chooseAlgo<<endl;
    cout<<"candidates file: "<<icf<<endl;
    cout<<"res file: "<<of<<endl;
    //advDP_* is index free
    if(4 != chooseAlgo || 3 != chooseAlgo){
        readIndex(infileIDX);
    }
	readGraph(infileGR);
    outOpenFile = fopen(outfile, "w");
    // run algos
    // BnB: branch and bound
    if(1 == chooseAlgo){
        run_BnB(srcNsink_num, queryFile, budget);
    }
    // Basic DP 
    if(2 == chooseAlgo){
        run_BasicDP(srcNsink_num, queryFile, budget);
    }
    // Adv DP is index free
    if(3 == chooseAlgo){
        run_AdvDP(srcNsink_num, queryFile, budget);
    }
    // Run Adv DP given candidates
    if(4 == chooseAlgo){ 
        run_AdvDP_given_candidate(srcNsink_num, queryFile, incandfile, budget);
    }
    // App and SimApp
    if(5 == chooseAlgo){
        run_appx(srcNsink_num, queryFile, icf, budget);
    }
    // bet/mst/pcp
    if(6 == chooseAlgo){
        run_apx_cand(srcNsink_num, queryFile, icf, budget);    
    }
    fclose(outOpenFile);
    return 0;
}