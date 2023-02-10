#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <cstring>
#include <cmath>
#define N 3
#define eps 0.00001

using namespace Eigen;
using namespace std;

int strcount (string str) {
	int n = -1;
        char c;
        ifstream fin(str);
        while (!fin.eof()) {
                fin.get(c);
                if (c == '\n')
                        n++;
        }
	fin.close();
	return n;
}


template <typename T>
vector< vector<T>> inarr (int n, int m, string str) {
        ifstream fin(str);
	vector <vector <T>> a(n, vector <T> (m));
	for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                        fin >> a[i][j];
        fin.close();
	return a;
}

template <typename T>
vector <T> inarr1 (int n, string str) {
        ifstream fin(str);
        vector <T> a(n);
        for (int i = 0; i < n; i++)
                        fin >> a[i];
        fin.close();
        return a;
}

template <typename S>
S maxInColumn (vector <vector <S>> a, int n, int col) {
	S maxC = a[0][col];
	for (int i = 0; i < n; i++) 
		if (a[i][col] > maxC)
			maxC = a[i][col];
	return maxC;
}

template <typename S>
S minInColumn (vector <vector <S>> a, int n, int col) {
        S minC = a[0][col];
        for (int i = 0; i < n; i++)
                if (a[i][col] < minC)
                        minC = a[i][col];
        return minC;
}


template <typename R>
vector <R> sameRemove (vector <R> a) {
        vector <R> b;
        b.push_back(a[0]);
        int t = 0;
        for (int i = 1; i < a.size(); i++) {
                for (int j = 0; j < b.size(); j++)
                        if (a[i] != b[j]) t++;
                if (t == b.size())
                        b.push_back(a[i]);
                t = 0;
        }
        return b;
}

vector <int> xBorderSearch (int n, int NE, int NP, vector <vector <int>> e, vector <vector <double>> nK, vector <int> nN) {
	vector <int> a;
	double minC = minInColumn <double> (nK, NP, n);
	double maxC = maxInColumn <double> (nK, NP, n);
	for (int i = 0; i < NP; i++)
		if (abs(nK[i][n] - maxC) < eps or abs(nK[i][n] - minC) < eps)
			for (int j = 0; j < NE; j++)
				for (int k = 0; k < 3; k++)
					if (e[j][k] == nN[i])
						a.push_back(j);
		

	vector <int> b = sameRemove <int> (a);
	return b;	
}

int main() {
/*
	string str_el;
	string str_nodN;
	string str_nodK;
	getline(cin, str_el);
	getline(cin, str_nodN);
	getline(cin, str_nodK);
*/
	double F = 1;
//	cin >> F;

	int NE = strcount("elements1.txt");
	int NP = strcount("nodesN.txt");
 	vector <vector <int>> elements(NE, vector <int> (N));
	vector <int> nodesN(NP);
        vector <vector <double>> nodesK(NP, vector <double> (2));

	elements = inarr <int> (NE, N, "elements1.txt");
	nodesN = inarr1 <int> (NP, "nodesN.txt");
	nodesK = inarr <double> (NP, 2, "nodesK.txt");
	vector <int> xBorder = xBorderSearch (0, NE, NP, elements, nodesK, nodesN);	

	return 0;
}




`






