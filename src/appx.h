#pragma once
#include "MSPAUtilities.h"
#include "H2H.h"
#include "pruned_landmark_labeling.h"

long long prim(const vector<vector<pair<int, int>>>& graph, int n) {
    // int n = graph.size();
    // vector<int> key(n, INF);
    vector<int> key(n, 2147483647);
    vector<bool> inMST(n, false);
    
    int totalWeight = 0;
    key[0] = 0;
    
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, 0});
    
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        
        if (inMST[u]) continue;
        inMST[u] = true;
        
        for (const auto& edge : graph[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (!inMST[v] && weight < key[v]) {
                key[v] = weight;
                pq.push({key[v], v});
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        totalWeight += key[i];
    }
    
    return totalWeight;


}
void fillG(vector<vector<pair<int, int>>>& newG, vector<int>& rns){
    newG.resize(rns.size());
    for(int i=0; i<rns.size(); i++){
        for(int j=0; j<rns.size(); j++){
            newG[i].push_back(make_pair(j, distanceQuery(rns[i], rns[j])));
        }
    }
}
void genG(vector<vector<pair<int, int>>>& newG, \
vector<pair<int, pair<vector<int>, long long>>> term_assndSrc){
    newG.resize(term_assndSrc.size());
    vector<int> vi;
    for(auto ele:term_assndSrc){
        vi.push_back(ele.first);
    }
    fillG(newG, vi);
}
//@deprecated replaced with index facilitated algo in gen bet cands.
void dij( vector<int>& termDist, vector<vector<int>>& pre, vector<int>& npath, int src){
    // init distances
    npath[src] = 1;
    termDist[src] = 0;
    auto initPath = [&](int u, int v) {
            pre[v] = {u};
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
                pre[nei].push_back(u);
                npath[nei] += npath[u];
            }
        }
    }

}
long long terminals_cand_tot_dis(vector<int>& terminal, int cand){
    long long t_tot_dis = 0;
    for(auto ele:terminal){
        t_tot_dis += distanceQuery(ele, cand);
    }
    return t_tot_dis;
}
long long assign_src2cand(vector<int>& terminal, vector<int>& aps, \
vector<vector<pair<int, int>>>& newG, unordered_map<int, vector<int>>& umpii){
    long long tot_dis = 0;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    for(int i=0; i<src.size(); i++){
        int prevDis = INT_MAX;
        int src2t = distanceQuery(src[i], sink);
        for(int j=0; j<aps.size(); j++){
            int tmpsrc2ap = distanceQuery(src[i], aps[j]);
            if(prevDis > tmpsrc2ap){
                prevDis = tmpsrc2ap;
            }
            if(prevDis > src2t){
                tot_dis += src2t;
            }else{
                tot_dis += prevDis;
            }
        }
    }
    aps.push_back(sink);
    fillG(newG, aps);
    aps.pop_back();
    return tot_dis;
}
// split only no assignment
// comp the tot tree weight
long long choose_root_cand(int root, int cand, vector<int>& src2split){
    long long choose_tot = distanceQuery(root, cand);
    for(int i=0; i<src2split.size(); i++){
        int toRoot = distanceQuery(src2split[i], root);
        int toCand = distanceQuery(src2split[i], cand);
        if(toRoot > toCand){
            choose_tot += toCand;
        }else{
            choose_tot += toRoot;
        }
    }
    return choose_tot;
}
//dfs grdily chooses the best and descent.
int find_cand_dfs(int root, long long& tot, vector<int>& src2split){
    int min_cand = -1;
    // Function to perform Depth-First Search on a graph
    vector<bool> visited(GraphSize+1, 0);
    stack<int> s;
    // Start DFS from the given vertex
    s.push(root);
    long long minDist = tot;
    int cnt=0;
    int decCnt=0;
    while (!s.empty()) {
        int curVer = s.top();
        s.pop();
        if(visited[curVer]) continue;
        visited[curVer] = 1;
        int bestNei = -1;
        long long bestTot = minDist;
        // Explore neighbors of the current vertex
        for (int i=0; i<Degree[curVer]; i++) {
            int nei = Neighbor[curVer][i];
            //@update 23.12.04 do not recalc visited ver
            if (visited[nei]) continue;
            long long tmpTot = choose_root_cand(root, nei, src2split); 
            cnt++;
            if(bestTot > tmpTot){
                bestTot = tmpTot;
                bestNei = nei;
            }
        }
        if(minDist > bestTot){
            decCnt++;
            minDist = bestTot;
            min_cand = bestNei;
        }else{
            break;
        }
        if(-1 == bestNei) continue;
        s.push(bestNei);
    }
    cout<<"dfs tot vis: "<<cnt<<endl;
    cout<<"dfs dec cnt: "<<decCnt<<endl;
    cout<<"dfs tot decr diff: "<<(tot-minDist)<<endl;
    tot = minDist;
    return min_cand;
}
//@udpated 23.12.05
//The immediate neis of a given root in bfs may not be satifying.
//We therefore init bfs/dfs with a cand.
// int init_fs(int root, long long& tot, vector<int>& src2split, vector<int>&exclude){
int init_fs(int root, long long& tot, vector<int>& src2split){
    int min_cand = -1;
    long long upBd = tot;
    long long minDist = 0;
    for (auto ele : src2split)
    {
        minDist += distanceQuery(ele, root);
    }
    vector<int> vis(GraphSize+1, 0);
    int cnt = 0;
    for (int i = 1; i < GraphSize + 1; i++)
    {
        if (i == root)
            continue;
        if (vis[i]) continue;
        long long tmpTot = choose_root_cand(root, i, src2split);
        // if (tmpTot > upBd){
        //     exclude[i] = 1;
        // }
        if (minDist >  tmpTot)
        {
            minDist = tmpTot;
            min_cand = i;
            cnt++;
            if(cnt > src2split.size()/2){
                break;
            }
        }
    } 
    tot = minDist;
    return min_cand;
}
//We use a smaller split tree to replace a subtree of prev root.
//bfs allow search on every deccresasing ver.
int find_cand_bfs(int root, long long& tot, vector<int>& src2split){
    int min_cand = -1;
    priority_queue<pair<long long, int>> pq;
    long long minDist = tot;
    int cnt = 0;
    int dffCnt = 0;
    long long diff = 0;
    vector<int> vis(GraphSize+1, 0);
    //speed up on social networks
    // vector<int> enQed(GraphSize+1, 0);

    pq.emplace(make_pair((-1*tot), root));
    vis[root] = 1;
    // enQed[root] = 1;
    while(!pq.empty()){
        int u = pq.top().second;
        long long uTot = pq.top().first*-1;
        pq.pop();
        // if(1 == vis[u]) continue;
        vis[u] = 1;
        long long tmpDiff = tot-uTot;
        if(tmpDiff > diff){
            diff = tmpDiff;
            minDist = uTot;
            min_cand = u;
            dffCnt++;
        }
        else{
            //kind of hybrid of bfs and dfs, we ignore unpromising ver's neighbors to further speedup.
            continue;
        }
        for(int i=0; i<Degree[u]; i++){
            int nei = Neighbor[u][i];
            if(1 == vis[nei]) continue;
            // if(1 == enQed[nei]) continue;
            long long tmpTot = choose_root_cand(root, nei, src2split);
            //prune unpromising nodes
            if(tmpTot >= tot) continue;
            pq.emplace(make_pair((-1*tmpTot), nei));
            vis[nei] = 1;
            // enQed[nei] = 1;
            cnt++;
        }
    }
    cout<<"tot vis "<<cnt<<endl;
    cout<<"tot dec: "<<dffCnt<<endl;
    cout<<"max diff: "<<diff<<endl;
    tot = minDist;
    return min_cand;
}
void allocate_cand_root(int u, int min_cand, vector<int>& src2split, \
vector<int>& split2cand, vector<int>& split2root, long long& cand_tot, long long& root_tot){
    for(int i=0; i<src2split.size(); i++){
        int toRoot = distanceQuery(src2split[i], u);
        int toCand = distanceQuery(src2split[i], min_cand);
        if(toRoot > toCand){
            cand_tot += toCand;
            split2cand.push_back(src2split[i]);
        }else{
            root_tot += toRoot;
            split2root.push_back(src2split[i]);
        }
    }
}

// With the help of indices, we enumerate nodes on the shortest paths from src to sink.
// We cnt the times each nodes is viewed and sort by cnt.
// The most viewed vers are the top ranked bet cands.
void gen_bet_with_idx(vector<int>& terminal, vector<int>& cand){
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<pair<int, int>> betCnt;
    for(int i=0; i<=GraphSize; i++){
        betCnt.push_back(make_pair(i, 0));
    }
    for(auto ele:srcs){
        int guideDis = distanceQuery(ele, sink);
        vector<int> vis(GraphSize+1, 0);
        queue<int> qi;
        qi.emplace(sink);
        while(!qi.empty()){
            int u = qi.front();
            qi.pop();
            if(vis[u]) continue;
            vis[u] = 1;
            for(int kth=0; kth<Degree[u]; kth++){
                int nei = Neighbor[u][kth];
                // if(vis[nei]) continue;
                int src2nei = distanceQuery(sink, nei); 
                // int src2nei = Weight[u][kth];
                //we start from sink.
                int nei2sink = distanceQuery(nei, ele);
                // int nei2sink = distanceQuery(sink, nei);
                if(guideDis == src2nei+nei2sink){
                    qi.push(nei);
                    betCnt[nei].second += 1;
                } 
            }
        }
    }
    betCnt[sink].second = 1;
    sort(betCnt.begin(), betCnt.end(), [&](pair<int, int> n1, pair<int, int> n2){
        return n1.second > n2.second; 
    });
    for(int i=0; i<10; i++){
        cout<<"["<<betCnt[i].first<<", "<<betCnt[i].second<<"]";
    }
    cout<<endl;
    int cnt=0;
    // k<|S|
    while(cand.size() < terminal.size()){
        if(betCnt[cnt].second == 0) break;
        cand.push_back(betCnt[cnt].first);
        cnt++;
    }
}

//prim based
//Kind of like random guess on graph.
//Every ver of G is on mst, therefore a
void gen_mst(vector<int>& terminal, vector<int>& cand){
    int n = GraphSize+1;
    vector<int> key(n, INT_MAX);
    vector<bool> inMST(n, false);
    //inpos starts from 1
    vector<int> inPos(n, 0);
    vector<pair<int, long long>> inTree;

    int srcCnt = 0;
    unordered_set<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    priority_queue<pair<int, int>, std::vector<pair<int, int>>, std::greater<pair<int, int>>> pq;
    //start from sink as the root of mst.
    pq.push({0, sink});
    key[sink]=0;
    inTree.push_back(make_pair(sink, terminals_cand_tot_dis(terminal, sink)));
    inPos[sink] = inTree.size();
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if(inMST[u]) continue;
        inMST[u] = true;
        // Include the popped vertex in MST
        if(srcs.find(u) != srcs.end()){
            srcCnt++;
        }
        if(srcCnt>=srcs.size()){
            break;
        }
        // Traverse all adjacent vertices of the popped vertex
        for (int i=0; i<Degree[u]; i++) {
            int nei = Neighbor[u][i];
            int wei = Weight[u][i];
            // Update the key and parent if v is not in MST and the weight is smaller
            if (!inMST[nei] && wei < key[nei]) {
                key[nei] = wei;
                pq.push({key[nei], nei});
                // cand.push_back(nei);
                if(inPos[nei] == 0){
                    inTree.push_back(make_pair(nei, terminals_cand_tot_dis(terminal, nei)));
                    inPos[nei] = inTree.size();
                }else{
                    inTree[inPos[nei]-1] = make_pair(nei, terminals_cand_tot_dis(terminal, nei));
                }
            }
        }
    }
    sort(inTree.begin(), inTree.end(), [&](pair<int, long long> pa, pair<int, long long> pb){
        return pa.second<pb.second;
    });
    for(int i=0; i<10; i++){
        cout<<"["<<inTree[i].first<<", "<<inTree[i].second<<"]";
    }
    cout<<endl;
    // k<|S|
    for(int i=0; i<terminal.size(); i++){
        cand.push_back(inTree[i].first);
    }
}

void gen_pcp_with_idx(vector<int>& terminal, vector<int>& cand){
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<pair<int, long long>> inPath;
    vector<int> inPos(GraphSize+1, 0);
    for(auto ele:srcs){
        int guideDis = distanceQuery(ele, sink);
        vector<int> vis(GraphSize+1, 0);
        queue<int> qi;
        qi.push(ele);
        inPath.push_back(make_pair(ele, terminals_cand_tot_dis(terminal, ele)));
        inPos[ele] = inPath.size();
        while(!qi.empty()){
            int u = qi.front();
            qi.pop();
            if(vis[u]) continue;
            vis[u] = 1;
            for(int kth=0; kth<Degree[u]; kth++){
                int nei = Neighbor[u][kth];
                int src2nei = distanceQuery(ele, nei);
                int nei2sink = distanceQuery(nei, sink);
                if(guideDis == src2nei+nei2sink){
                    qi.push(nei);
                    if(inPos[nei]  == 0){
                        inPath.push_back(make_pair(nei, terminals_cand_tot_dis(terminal, nei)));
                        inPos[nei] = inPath.size();
                    }else{
                        inPath[inPos[nei]-1] = make_pair(nei, terminals_cand_tot_dis(terminal, nei));
                    }
                    break;//break for loop, as only one sp is needed
                } 
            }
        }
    }
    // valued by avg dis to
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // shuffle(tmp.begin(), tmp.end(), gen);
    sort(inPath.begin(), inPath.end(), [&](pair<int, long long> pa, pair<int, long long> pb){
        return pa.second<pb.second;
    });
    for(int i=0; i<10; i++){
        cout<<"["<<inPath[i].first<<", "<<inPath[i].second<<"]";
    }
    cout<<endl;
    // k<|S|
    for(int i=0; i<terminal.size(); i++){
        cand.push_back(inPath[i].first);
    }
}