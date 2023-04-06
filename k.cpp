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

int maxInArr (vector <int> arr) {
	int t = arr[0];
	for (int i = 0; i < arr.size(); i++)
		if (t < arr[i])
			t = arr[i];
	return t;
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

vector <int> xBorderSearch (int n, int NE, int NP, vector <vector <int>> e, vector <vector <float>> nK, vector <int> nN) {
	vector <int> a;
	float minC = minInColumn <float> (nK, NP, n);
	float maxC = maxInColumn <float> (nK, NP, n);
	for (int i = 0; i < NP; i++)
		if (abs(nK[i][n] - maxC) < eps or abs(nK[i][n] - minC) < eps)
			for (int j = 0; j < NE; j++)
				for (int k = 0; k < 3; k++)
					if (e[j][k] == nN[i])
						a.push_back(j);
		

	vector <int> b = sameRemove <int> (a);
	return b;	
}

int nodeSearch (int m, vector <int> arr) {
	int i;
	for (i = 0; i < arr.size(); i++) 
		if (arr[i] == m)
			break;
	return i;

}

class Element {
public:	
	Vector3f x;
	Vector3f y;
	Vector3f nod;
	int f;
	Matrix3f C;
	Matrix <float, 3, 6> B;
	Matrix <float, 6, 6> k;
};

int main() {
/*
	string str_el;
	string str_nodN;
	string str_nodC;
	getline(cin, str_el);
	getline(cin, str_nodN);
	getline(cin, str_nodC);
*/
	float F = 1;
	float Young = 1200;
	float Poisson = 0.2;
//	cin >> F;
//	cin >> Young;
//	cin >> Poisson

	int NE = strcount("elements1.txt");
	int NP = strcount("nodesN.txt");
 	vector <vector <int>> elements(NE, vector <int> (N));
	vector <int> nodesN(NP);
        vector <vector <float>> nodesC(NP, vector <float> (2));

	elements = inarr <int> (NE, N, "elements1.txt"); //elements numbers NExN
	nodesN = inarr1 <int> (NP, "nodesN.txt"); //nodes numbers NP
	nodesC = inarr <float> (NP, 2, "nodesC.txt"); //coordinates of nodes NPx2
	vector <int> xBorder = xBorderSearch (0, NE, NP, elements, nodesC, nodesN);

	Matrix3f D;
	D <<   	1, Poisson, 0,
    		Poisson, 1, 0,
    		0, 0, (1 - Poisson) / 2;

	D *= Young / (1 - pow(Poisson, 2));

        SparseMatrix <float> K(2*maxInArr(nodesN) + 2, 2*maxInArr(nodesN) + 2);
        vector <Triplet <float>> tripl;


	vector <Element> element(NE);
	for (int i = 0; i < NE; i++) {
		element[i].nod << elements[i][0], elements[i][1], elements[i][2];
		vector <float> temp(6);
		for (int j = 0; j < 3; j++) {
			int t = nodeSearch(elements[i][j], nodesN);
			temp[2*j] = nodesC[t][0];
			temp[2*j + 1] = nodesC[t][1];
		}

		element[i].x << temp[0], temp[2], temp[4];
		element[i].y << temp[1], temp[3], temp[5];
		element[i].C << Vector3f(1, 1, 1), element[i].x, element[i].y;
		Matrix3f iC = element[i].C.inverse();
		for (int j = 0; j < 3; j++) {

			element[i].B(0, 2*j) = iC(1, j);
			element[i].B(0, 2*j + 1) = 0;
			element[i].B(1, 2*j) = 0;
			element[i].B(1, 2*j + 1) = iC(2, j);
			element[i].B(2, 2*j) = iC(2, j);
			element[i].B(2, 2*j + 1) = iC(1, j);
		}
		element[i].k = element[i].B.transpose() * D * element[i].B * element[i].C.determinant() / 2;
		
		for (int j = 0; j < 3; j++)
			for (int l = 0; l < 3; l++) {

			Triplet <float> trpl11(2*element[i].nod[j], 2*element[i].nod[l], element[i].k(2*j, 2*l));
			Triplet <float> trpl12(2*element[i].nod[j], 2*element[i].nod[l] + 1, element[i].k(2*j, 2*l + 1));
			Triplet <float> trpl21(2*element[i].nod[j] + 1, 2*element[i].nod[l], element[i].k(2*j + 1, 2*l));
			Triplet <float> trpl22(2*element[i].nod[j] + 1, 2*element[i].nod[l] + 1, element[i].k(2*j + 1, 2*l + 1));

			tripl.push_back(trpl11);
			tripl.push_back(trpl12);
			tripl.push_back(trpl21);
			tripl.push_back(trpl22);
		}
	
	}

	K.setFromTriplets(tripl.begin(), tripl.end());
	K.makeCompressed();
	cout << K;

	return 0;
}

