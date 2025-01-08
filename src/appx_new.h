#pragma once
#include <algorithm>
#include "MSPAUtilities.h"
#include "H2H.h"
#include "pruned_landmark_labeling.h"

long long prim(const vector<vector<pair<int, int>>>& graph, int n) {
    // int n = graph.size();
    vector<int> key(n, INF);
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
void fillG(vector<vector<pair<int, int>>>& newG, vector<int>& rns, int idx_flag){
    if(idx_flag){
        fillG(newG, rns);
    }else{
        newG.resize(rns.size());
        for(int i=0; i<rns.size(); i++){
            for(int j=0; j<rns.size(); j++){
                newG[i].push_back(make_pair(j, dij(rns[i], rns[j])));
            }
        }

    }
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

long long assign_src2cand(vector<int>& terminal, vector<int>& aps, \
vector<vector<pair<int, int>>>& newG, unordered_map<int, vector<int>>& umpii, int idx_flag){
    long long tot_dis = 0;
    int sink = terminal.back();
    vector<int> src(terminal.begin(), terminal.end()-1);
    for(int i=0; i<src.size(); i++){
        int prevDis = INT_MAX;
        
        int src2t;
        if(idx_flag){
            src2t = distanceQuery(src[i], sink);
        }else{
            src2t = dij(src[i], sink);
        }
        for(int j=0; j<aps.size(); j++){
            int tmpsrc2ap;
            if(idx_flag){
                tmpsrc2ap = distanceQuery(src[i], aps[j]);
            }else{
                tmpsrc2ap = dij(src[i], aps[j]);
            }
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
    fillG(newG, aps, idx_flag);
    aps.pop_back();
    return tot_dis;
}
long long choose_root_cand(int root, int cand, vector<int>& src2split, int idx_flag){
    long long choose_tot = 0;
    for(int i=0; i<src2split.size(); i++){
        int toRoot;
        int toCand;
        if(idx_flag){
            toRoot = distanceQuery(src2split[i], root);
            toCand = distanceQuery(src2split[i], cand);
        }else{
            toRoot = dij(src2split[i], root);
            toCand = dij(src2split[i], cand);
        }
        if(toRoot > toCand){
            choose_tot += toCand;
        }else{
            choose_tot += toRoot;
        }
    }
    return choose_tot;
}
int tmp_cand_dfs(int u, long long& tot, vector<int>& src2split){
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
            long long tmpTot = choose_root_cand(root, nei, src2split, 1); 
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

int find_cand_dfs(int u, long long& tot, vector<int>& src2split){
    int min_cand = -1;
    // Function to perform Depth-First Search on a graph
    int start=u;
    vector<bool> visited(GraphSize+1, 0);
    stack<int> s;
    // Start DFS from the given vertex
    s.push(start);
    visited[start] = 1;
    long long minDist = LLONG_MAX;
    int cnt=0;
    while (!s.empty()) {
        int currentVertex = s.top();
        s.pop();
        int bestNei = 0;
        long long bestTot = LLONG_MAX;
        // Explore neighbors of the current vertex
        for (int i=0; i<Degree[currentVertex]; i++) {
            int nei = Neighbor[currentVertex][i];
            long long tmpTot = choose_root_cand(currentVertex, nei, src2split, 1); 
            cnt++;
            if(bestTot > tmpTot){
                bestTot = tmpTot;
                bestNei = nei;
            }
        }
        if(minDist > bestTot){
            minDist = bestTot;
            min_cand = bestNei;
        }else{
            break;
        }
        s.push(bestNei);
    }
    cout<<"dfs tot vis: "<<cnt<<endl;
    tot = minDist;
    return min_cand;
}
int find_cand_bfs(int root, long long& tot, vector<int>& src2split){
    int min_cand = -1;
    priority_queue<pair<long long, int>> pq;
    long long minDist = tot;
    int cnt = 0;
    int dffCnt = 0;
    long long diff = 0;
    vector<int> vis(GraphSize+1, 0);
    vector<int> enQed(GraphSize+1, 0);

    pq.emplace(make_pair((-1*tot), root));
    enQed[root] = 1;
    while(!pq.empty()){
        int u = pq.top().second;
        long long uTot = pq.top().first*-1;
        pq.pop();
        // if(uTot > minDist){
            // break;
        // }else{
            // minDist = uTot;
            // min_cand = u;
        // }
        if(1 == vis[u]) continue;
        vis[u] = 1;
        long long tmpDiff = tot-uTot;
        if(tmpDiff > diff){
            diff = tmpDiff;
            minDist = uTot;
            min_cand = u;
            dffCnt++;
        }

        for(int i=0; i<Degree[u]; i++){
            int nei = Neighbor[u][i];
            if(1 == vis[nei]) continue;
            if(1 == enQed[nei]) continue;
            long long tmpTot = choose_root_cand(u, nei, src2split, 1);
            if(tmpTot > tot) continue;
            pq.emplace(make_pair((-1*tmpTot), nei));
            enQed[nei] = 1;
            cnt++;
        }
    }
    cout<<"tot vis "<<cnt<<endl;
    cout<<"tot dec: "<<dffCnt<<endl;
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
void gen_bet(vector<int>& terminal, vector<int>& cand){
    // void dij( vector<int>& termDist, vector<vector<int>>& pre, vector<int>& npath, int src){
    vector<int> nodes;
    for(int i=1; i<=GraphSize+1; i++){
        nodes.push_back(i);
    }
    vector<int> termDist(GraphSize+1, INF); 
    vector<vector<int>> pre; 
    pre.resize(0);
    vector<int> npath(GraphSize+1, 0);
    //collect number of shortest paths pass through by from each srcs
    vector<int> tmpNPath;
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    for(auto ele:srcs){
        // dij(termDist, pre, tmpNPath, ele);
        dij(termDist, pre, tmpNPath, ele);
        for(int i=0; i<GraphSize+1; i++){
            npath[i]+=tmpNPath[i];
        }
        vector<int> (GraphSize+1, INF).swap(termDist);
        vector<int> (GraphSize+1, 0).swap(tmpNPath);
    }
    //gt for descending
    sort(nodes.begin(), nodes.end(), [&](int n1, int n2){
        return npath[n1]>npath[n2];
    });
    int cnt=0;
    while(cand.size() < 10000){
        if(nodes[cnt] == 0 ) break;
        cand.push_back(nodes[cnt]);
        cnt++;
    }
}

void gen_bet_with_idx(vector<int>& terminal, vector<int>& cand){
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<int> betCnt(GraphSize+1, 0); 
    vector<int> nodes;
    for(int i=1; i<=betCnt.size(); i++){
        nodes.push_back(i);
    }
    for(auto ele:srcs){
        int guideDis = distanceQuery(ele, sink);
        vector<int> vis(GraphSize+1, 0);
        queue<int> qi;
        qi.emplace(sink);
        betCnt[ele] += 1;
        betCnt[sink]+= 1;
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
                    betCnt[nei]++;
                } 
            }
        }
    }
    //most viewed is the 1st elem
    //gt for descending
    sort(nodes.begin(), nodes.end(), [&](int n1, int n2){
        return betCnt[n1]>betCnt[n2];
    });
    int cnt=0;
    while(cand.size() < 10000){
        if(nodes[cnt] == 0 ) break;
        cand.push_back(nodes[cnt]);
        cnt++;
    }
}

// Prim's algorithm to find MST
void gen_mst(vector<int>& terminal, vector<int>& cand) {
    int n = GraphSize+1;
    // Create a vector to store the key (weight) of each vertex
    std::vector<int> key(n, INT_MAX);
    // Create a vector to mark whether a vertex is included in MST
    std::vector<bool> inMST(n, false);

    int srcCnt = 0;
    unordered_set<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    std::priority_queue<pair<int, int>, std::vector<pair<int, int>>, std::greater<pair<int, int>>> pq;
    pq.push({0, sink});
    key[sink]=0;
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
                cand.push_back(nei);
            }
        }
    }
    // gt for descending
    sort(cand.begin(), cand.end(), [&](int v1, int v2){
        return Degree[v1] > Degree[v2];
    });
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // shuffle(cand.begin(), cand.end(), gen);
}

void gen_pcp(vector<int>& terminal, vector<int>& cand){
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<int> tmp;
    set<int> tmpsi;
    vector<int> termDist(GraphSize+1, INF); 
    vector<int> npath(GraphSize+1, 0);
    vector<vector<int>> pre;
    pre.resize(GraphSize+1);
    for(auto ele:srcs){
        dij(termDist, pre, npath, ele);
        int cur = sink;
        while(cur != ele){
            if(tmpsi.find(cur) == tmpsi.end()){
                tmp.push_back(cur);
            }
            cur=pre[cur].back();
        } 
        vector<int>(GraphSize+1, INF).swap(termDist);
        vector<vector<int>> (GraphSize+1, vector<int>()).swap(pre);
        vector<int> (GraphSize+1, 0).swap(npath);
    }

    //gt for descending
    sort(tmp.begin(), tmp.end(), [&](int v1, int v2){
        return Degree[v1] > Degree[v2];
    });
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // shuffle(tmp.begin(), tmp.end(), gen);
    // int candNum; 
    // if(tmp.size() > 10000){
    //     candNum = 10000;
    // }else{
    //     candNum = tmp.size();
    // }
    for(int i=0; i<tmp.size(); i++){
        cand.push_back(tmp[i]);
    } 
}

void gen_pcp_with_idx(vector<int>& terminal, vector<int>& cand){
    vector<int> srcs(terminal.begin(), terminal.end()-1);
    int sink = terminal.back();
    vector<int> tmp;
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
                int src2nei = distanceQuery(ele, nei);
                int nei2sink = distanceQuery(nei, sink);
                if(guideDis == src2nei+nei2sink){
                    qi.push(nei);
                    tmp.push_back(nei);
                    break;
                } 
            }
        }
    }
    //gt for descending
    sort(tmp.begin(), tmp.end(), [&](int v1, int v2){
        return Degree[v1] > Degree[v2];
    });
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // shuffle(tmp.begin(), tmp.end(), gen);
    // int candNum; 
    // if(tmp.size() > 10000){
    //     candNum = 10000;
    // }else{
    //     candNum = tmp.size();
    // }
    for(int i=0; i<tmp.size(); i++){
        cand.push_back(tmp[i]);
    }
}