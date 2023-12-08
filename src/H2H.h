#pragma once
#include <algorithm>
#include <bitset>
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
#include <sys/time.h>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_set>
#include <xmmintrin.h>
#include "pruned_landmark_labeling.h"
using namespace std;
const int infinity = 999999999;
const int SIZEOFINT = 4;

PrunedLandmarkLabeling<64> pll;
int *toRMQ, *height, **RMQIndex;
int *belong;
int root, TreeSize;
int **rootToRoot, *rootSite;
int **dis, **pos, **pos2;
int *posSize, *pos2Size;
int *chSize;
int **ch;
int *LOG2, *LOGD;

int *Degree;
int **Neighbor, **Weight;
int GraphSize;
int EdgeSize;

int n;
int *EulerSeq;


FILE* fin;

double GetTime(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}
void scanIntArray(int *a, int n)
{
	fread(a, SIZEOFINT, n, fin);
}
int *scanIntVector(int *a)
{
	int _n;
	fread(&_n, SIZEOFINT, 1, fin);
	a = (int *)malloc(sizeof(int) * _n);
	scanIntArray(a, _n);
	return a;
}
//get adjList of G
//in Neighbor[][]; Weight[][]; Degree[]
void readGraph(char *filename)
{
		FILE *file = fopen(filename, "r");
		int n, m;
		fscanf(file, "%d %d", &n, &m);
		GraphSize = n;
		EdgeSize = m;
		cout << "GraphSize: " << GraphSize << endl;
		cout << "EdgeSize: " << EdgeSize << endl;
		Degree = (int *)malloc(sizeof(int) * (n + 1));
		vector<vector<pair<int, int>>> nb;
		vector<pair<int, int>> v;
		v.clear();
		for (int i = 0; i <= n; i++)
		{
			nb.push_back(v);
		}
		for (int i = 0; i < m; i++)
		{
			int x, y, z;
			fscanf(file, "%d %d %d", &x, &y, &z);
			nb[x].push_back(make_pair(y, z));
		}
		Neighbor = (int **)malloc(sizeof(int *) * (n + 1));
		Weight = (int **)malloc(sizeof(int *) * (n + 1));
		for (int i = 1; i <= n; i++)
		{
			Degree[i] = nb[i].size();
			Neighbor[i] = (int *)malloc(sizeof(int) * nb[i].size());
			Weight[i] = (int *)malloc(sizeof(int) * nb[i].size());
			for (int j = 0; j < (int)nb[i].size(); j++)
			{
				Neighbor[i][j] = nb[i][j].first;
				Weight[i][j] = nb[i][j].second;
			}
		}
}

void readIndex(char *file)
{
	string ifstr = file;
	size_t strpos = ifstr.find("-d");
	if(strpos!=string::npos){
		double _time = GetTime();
		int tree_height = 0, tree_width = 0, most_sp = 0;
		fin = fopen(file, "rb");
		fread(&n, SIZEOFINT, 1, fin);
		int ts;
		fread(&ts, SIZEOFINT, 1, fin);
		TreeSize = ts;
		height = (int *)malloc(sizeof(int) * (ts + 1));
		for (int i = 0; i < ts; i++)
		{
			fread(&height[i], SIZEOFINT, 1, fin);
		}
		belong = (int *)malloc(sizeof(int) * (n + 1));
		fread(belong, SIZEOFINT, n + 1, fin);
		toRMQ = (int *)malloc(sizeof(int) * (n + 1));
		fread(toRMQ, SIZEOFINT, n + 1, fin);
		int ris;
		fread(&ris, SIZEOFINT, 1, fin);
		fread(&ts, SIZEOFINT, 1, fin);
		EulerSeq = (int *)malloc(sizeof(int) * (ts + 1));
		RMQIndex = (int **)malloc(sizeof(int *) * (ris + 1));
		for (int i = 0; i < ris; i++)
		{
			RMQIndex[i] = scanIntVector(RMQIndex[i]);
		}
		fread(&root, SIZEOFINT, 1, fin);
		cout << "root: " << root << endl;

		posSize = (int *)malloc(sizeof(int) * (n + 1));
		pos2Size = (int *)malloc(sizeof(int) * (n + 1));
		pos = (int **)malloc(sizeof(int *) * (TreeSize));
		pos2 = (int **)malloc(sizeof(int *) * (TreeSize));
		dis = (int **)malloc(sizeof(int *) * (TreeSize));
		chSize = (int *)malloc(sizeof(int) * (TreeSize));
		ch = (int **)malloc(sizeof(int *) * (TreeSize));

		for (int i = 0; i < TreeSize; i++)
		{
			fread(&chSize[i], SIZEOFINT, 1, fin);
			ch[i] = (int *)malloc(sizeof(int) * chSize[i]);
			for (int j = 0; j < chSize[i]; j++)
			{
				int x;
				fread(&x, SIZEOFINT, 1, fin);
				ch[i][j] = x;
			}
		}
		for (int i = 0; i < TreeSize; i++)
		{
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			fread(&posSize[x], SIZEOFINT, 1, fin);
			pos[x] = (int *)malloc(sizeof(int) * (posSize[x] + 1));
			fread(pos[x], SIZEOFINT, posSize[x], fin);
			if (posSize[x] > tree_width)
				tree_width = posSize[x];
			int _n;
			fread(&_n, SIZEOFINT, 1, fin);
			dis[x] = (int *)malloc(sizeof(int) * _n);
			fread(dis[x], SIZEOFINT, _n, fin);
			if (_n > tree_height)
				tree_height = _n;
		}
		printf("dis read finished!\n");
		for (int i = 0; i < TreeSize; i++)
		{
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			fread(&pos2Size[x], SIZEOFINT, 1, fin);
			pos2[x] = (int *)malloc(sizeof(int) * (pos2Size[x] + 1));
			fread(pos2[x], SIZEOFINT, pos2Size[x], fin);
			if (pos2Size[x] > most_sp)
				most_sp = pos2Size[x];
		}
		// initialize RMQ table____START
		LOG2 = (int *)malloc(sizeof(int) * (n * 2 + 10));
		LOGD = (int *)malloc(sizeof(int) * (n * 2 + 10));
		int k = 0, j = 1;
		for (int i = 0; i < n * 2 + 10; i++)
		{
			if (i > j * 2)
			{
				j *= 2;
				k++;
			}
			LOG2[i] = k;
			LOGD[i] = j;
		}
		// initialize RMQ table____END
		fclose(fin);
		printf("Load Index Time : %lf sec\n", (GetTime() - _time));
		printf("tree height: %d\n", tree_height);
		printf("tree width: %d\n", tree_width);
		printf("most search space: %d\n", most_sp);
	}else{
		pll.LoadIndex(file);
	}
}
inline int LCAQuery(int _p, int _q)
{
	int p = toRMQ[_p], q = toRMQ[_q];

	if (p > q)
	{
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;

	int i = LOGD[len], k = LOG2[len];

	q = q - i + 1;
	if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
		return RMQIndex[k][p];
	else
		return RMQIndex[k][q];
}
inline int distanceQuery(int p, int q)
{
	if(Weight[p][0] != 1){
		if (p == q)
			return 0;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y)
		{
			if (lca == y)
			{
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			return dis[y][pos[x][posSize[x] - 1]];
		}
		else
		{
			int res = infinity;
			int *dx = dis[x], *dy = dis[y], *p2 = pos2[lca];
			_mm_prefetch(dx, _MM_HINT_T0);
			_mm_prefetch(dy, _MM_HINT_T0);
			_mm_prefetch(p2, _MM_HINT_T0);
			int ps = pos2Size[lca];
			for (int i = 0; i < ps; i++)
			{
				int tmp = dx[p2[i]] + dy[p2[i]];
				if (res > tmp)
					res = tmp;
			}
			return res;
		}
	}
	else{
		// return pll.QueryDistance(p-1, q-1);
		return pll.QueryDistance(p, q);
	}
}