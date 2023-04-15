#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"
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

vector <int> BorderSearch (int n, int NE, int NP, vector <vector <int>> e, vector <vector <float>> nC, vector <int> nN) {
	vector <int> a;
	float minC = minInColumn <float> (nC, NP, n);
	float maxC = maxInColumn <float> (nC, NP, n);
	for (int i = 0; i < NP; i++) {
		int t = 0;
		if (abs(nC[i][n] - maxC) < eps or abs(nC[i][n] - minC) < eps)
//			a.push_back(nN[i]);
			for (int j = 0; j < NE; j++) {
				for (int k = 0; k < 3; k++)
					if (e[j][k] == nN[i])
						a.push_back(j);
			}
	}

	vector <int> b = sameRemove <int> (a);
	vector <int> c;
	for (int i = 0; i < b.size(); i++) {
		int t = 0;
		int u, v;
		for (int j = 0; j < 3; j++)
			if ((nC[e[b[i]][j]][n] == minC) or (nC[e[b[i]][j]][n] == maxC)) {
				t++;
				if (t == 1) u = e[b[i]][j];
				if (t == 2) v = e[b[i]][j];
			}
		if (t == 2) {
			c.push_back(u);
			c.push_back(v);
		}
	}

	return c;	
}

int maxInArr (vector <int> arr) {
        int t = arr[0];
        for (int i = 0; i < arr.size(); i++)
                if (t < arr[i])
                        t = arr[i];
        return t;
}

vector <float> IJCalc (float a, float b, float c, float d, int xi, int eta) {
	vector <float> ij;
	float det;
	det = b*c - a*d + b*c * xi - a*d * xi - b*c * eta + a*d * eta + b*c * xi*eta - a*d * xi*eta;
	ij.push_back(2*(c - d * xi*xi*xi + c * xi*xi + 2*d * xi*xi - 2*c * xi - d * xi) / det);
        ij.push_back(2*(- a + b * xi*xi*xi - a * xi*xi - 2*b * xi*xi + 2*a * xi + b * xi) / det);
        ij.push_back(2*(- c - d + c * xi + d * xi - c * eta + d * eta + c * xi*eta - d * xi*eta) / det);
        ij.push_back(2*(a + b - a * xi - b * xi + a * eta - b * eta - a * xi*eta + b * xi*eta) / det);
	return ij;
}



int nodeSearch (int m, vector <int> arr) {
	int i;
	for (i = 0; i < arr.size(); i++) 
		if (arr[i] == m)
			break;
	return i;

}

float Nxi (float eta, int i) {
	float Nn;
	if (i == 0)
		Nn = (-1 - eta) / 4;
	if (i == 1)
		Nn = (-1 + eta) / 4;
	if (i == 2) 
		Nn = (1 - eta) / 4;
	if (i == 3)
		Nn = (1 + eta) / 4;
	return Nn;
}

float Neta (float xi, int i) {
        float Nn;
	if (i == 0)
        	Nn = (1 - xi) / 4;
	if (i == 1)
		Nn = (-1 + xi) / 4;
	if (i ==  2)
		Nn = (-1 - xi) / 4;
	if (i == 3)
		Nn = (1 + xi) / 4;
        return Nn;
}



class Element {
public:	
	Vector3f x;
	Vector3f y;
	Vector3f nod;
	Matrix3f C;
	Vector <float, 6> u;
	Vector3f sigma;
	Vector3f epsilon;
	Matrix <float, 3, 6> B;
	Matrix <float, 6, 6> k;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class InfElement {
public:
	
	Vector <float, 4> P;
	Vector <float, 4> Q;

	Vector <float, 8> u;
	Vector3f sigma;
	Vector3f epsilon;
	Matrix <float, 3, 8> B;
	Matrix <float, 8, 8> k;
	
	vector <float> iJ; 
		
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
	float Young = 2000;
	float Poisson = 0.3;
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
	
//	vector <vector <int>> elements2 (NE, vector <int> (N));
//        vector <int> nodesN2 (NP);

	for (int i = 0; i < NP; i++) {
		for (int j = 0; j < NE; j++)
			for (int l = 0; l < 3; l++)
				if (elements[j][l] == nodesN[i])
					elements[j][l] = i;
		nodesN[i] = i;
	}	
 

	Matrix3f D;
	D <<   	1, Poisson, 0,
    		Poisson, 1, 0,
    		0, 0, (1 - Poisson) / 2;

	D *= Young / (1 - pow(Poisson, 2));
	
	SparseMatrix <float> K (2*NP, 2*NP);	
	vector <Triplet <float>> tripl;

	vector <int> xBorder = BorderSearch (0, NE, NP, elements, nodesC, nodesN);
        vector <int> yBorder = BorderSearch (1, NE, NP, elements, nodesC, nodesN);

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
// Global matrix		
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

// Infinite elements

	float Xc = (maxInColumn(nodesC, NP, 0) + minInColumn(nodesC, NP, 0)) / 2;
	float Yc = (maxInColumn(nodesC, NP, 1) + minInColumn(nodesC, NP, 1)) / 2;

	int NI = (xBorder.size() + yBorder.size()) / 2;
	
	vector <InfElement> infel (NI);
	
	for (int i = 0; i < xBorder.size() / 2; i++) {
		infel[i].P << nodesC[xBorder[2*i]][0], nodesC[xBorder[2*i]][1], nodesC[xBorder[2*i + 1]][0], nodesC[xBorder[2*i + 1]][1];
		infel[i].Q << 2*infel[i].P[0], 2*infel[i].P[1], 2*infel[i].P[2], 2*infel[i].P[3];
	}
	for (int i = 0; i < yBorder.size() / 2; i++) {
                infel[i + xBorder.size() / 2].P << nodesC[yBorder[2*i]][0], nodesC[yBorder[2*i]][1], nodesC[yBorder[2*i + 1]][0], nodesC[yBorder[2*i + 1]][1];
                infel[i + xBorder.size() / 2].Q << 2*infel[i + xBorder.size() / 2].P[0], 2*infel[i + xBorder.size() / 2].P[1],
						   2*infel[i + xBorder.size() / 2].P[2], 2*infel[i + xBorder.size() / 2].P[3];
	}
	

	vector <float> GPxi {-0.5774, -0.5774, 0.5774, 0.5774};
	vector <float> GPeta {0.5774, -0.5774, -0.5774, 0.5774};
	K.resize(2*NP + 2*NI, 2*NP + 2*NI);

	for (int i = 0; i < NI; i++) {
		
/*        	infel[i].J << 0.5 * ((eta / (1 - xi*xi)) * (infel[i].Q[0] - infel[i].Q[2]) + 1 / (1 - xi*xi) * (infel[i].Q[0] + infel[i].Q[2])),
              	     0.5 * (1 + xi / (1 - xi) * (infel[i].Q[0] - infel[i].Q[2])),
		     0.5 * ((eta / (1 - xi*xi)) * (infel[i].Q[1] - infel[i].Q[3]) + 1 / (1 - xi*xi) * (infel[i].Q[1] + infel[i].Q[3])),
                     0.5 * (1 + xi / (1 - xi) * (infel[i].Q[1] - infel[i].Q[3]));

		infel[i].det = infel[i].Q[1] * infel[i].Q[2] - infel[i].Q[0] * infel[i].Q[3] +
		               infel[i].Q[1] * infel[i].Q[2] * xi - infel[i].Q[0] * infel[i].Q[3] * xi + 
			       infel[i].Q[1] * infel[i].Q[2] * eta - infel[i].Q[0] * infel[i].Q[3] * eta + 
			       infel[i].Q[1] * infel[i].Q[2] * xi * eta - infel[i].Q[0] * infel[i].Q[3] *xi * eta;

		infel[i].IJ << 
		2*(infel[i].Q[1] - infel[i].Q[3] * xi*xi*xi + infel[i].Q[1] * xi*xi + 2*infel[i].Q[3] * xi*xi - 2*infel[i].Q[1] * xi - infel[i].Q[3] * xi) / infel[i].det,
		2*(- infel[i].Q[0] + infel[i].Q[2] * xi*xi*xi - infel[i].Q[0] * xi*xi - 2*infel[i].Q[2] * xi*xi + 2*infel[i].Q[0] * xi + infel[i].Q[2] * xi) / infel[i].det,
		2*(- infel[i].Q[1] - infel[i].Q[3] + infel[i].Q[1] * xi + infel[i].Q[3] * xi - infel[i].Q[1] * eta + infel[i].Q[3] * eta + infel[i].Q[1] * xi*eta - infel[i].Q[3] * xi*eta) / infel[i].det,
		2*(infel[i].Q[0] + infel[i].Q[2] - infel[i].Q[0] * xi - infel[i].Q[2] * xi + infel[i].Q[0] * eta - infel[i].Q[2] * eta - infel[i].Q[0] * xi*eta + infel[i].Q[2] * xi*eta) / infel[i].det;

*/	
		for (int j = 0; j < 4; j++) {
			if (j == 0) {
				vector <float> ij;
				ij = IJCalc(infel[i].Q[0], infel[i].Q[2], infel[i].Q[1], infel[i].Q[3], -0.5774, 0.5774);
				infel[i].iJ.push_back(ij[0]);
				infel[i].iJ.push_back(ij[1]);
				infel[i].iJ.push_back(ij[2]);
				infel[i].iJ.push_back(ij[3]);
			}
			 if (j == 1) {
				vector <float> ij;
                                ij = IJCalc(infel[i].Q[0], infel[i].Q[2], infel[i].Q[1], infel[i].Q[3], -0.5774, -0.5774);
                                infel[i].iJ.push_back(ij[0]);
                                infel[i].iJ.push_back(ij[1]);
                                infel[i].iJ.push_back(ij[2]);
                                infel[i].iJ.push_back(ij[3]);
                        }
			 if (j == 2) {
				vector <float> ij;
                                ij = IJCalc(infel[i].Q[0], infel[i].Q[2], infel[i].Q[1], infel[i].Q[3], 0.5774, -0.5774);
                                infel[i].iJ.push_back(ij[0]);
                                infel[i].iJ.push_back(ij[1]);
                                infel[i].iJ.push_back(ij[2]);
                                infel[i].iJ.push_back(ij[3]); 
                        }
			 if (j == 3) {
				vector <float> ij;
                                ij = IJCalc(infel[i].Q[0], infel[i].Q[2], infel[i].Q[1], infel[i].Q[3], 0.5774, 0.5774);
                                infel[i].iJ.push_back(ij[0]);
                                infel[i].iJ.push_back(ij[1]);
                                infel[i].iJ.push_back(ij[2]);
                                infel[i].iJ.push_back(ij[3]); 
                        }
		}

		for (int j = 0; j < 8; j++)
			for(int l = 0; l < 8; l++)
				infel[i].k(j, l) = 0;

		for (int j = 0; j < 4; j++) {
			for(int l = 0; l < 4; l++) {
				infel[i].B(0, 2*l) = Nxi(GPeta[j], l) * infel[i].iJ[4*j] + Neta(GPxi[j], l) * infel[i].iJ[4*j + 1];
				infel[i].B(0, 2*l + 1) = 0;
				infel[i].B(1, 2*l + 1) = Nxi(GPeta[j], l) * infel[i].iJ[4*j + 2] + Neta(GPxi[j], l) * infel[i].iJ[4*j + 3];
				infel[i].B(1, 2*l) = 0;
                		infel[i].B(2, 2*l) = Nxi(GPeta[j], l) * infel[i].iJ[4*j + 2] + Neta(GPxi[j], l) * infel[i].iJ[4*j + 3];
				infel[i].B(2, 2*l + 1) = Nxi(GPeta[j], l) * infel[i].iJ[4*j] + Neta(GPxi[j], l) * infel[i].iJ[4*j + 1];
			infel[i].k += infel[i].B.transpose() * D * infel[i].B;
			}
		
		}
	}

	for (int i = 0; i < NE + NI; i++) {
		if (i < NE)
		for (int j = 0; j < 3; j++)
                        for (int l = 0; l < 3; l++) {

                        Triplet <float> trpl11(2*element[i].nod[j], 2*element[i].nod[l], 0);
                        Triplet <float> trpl12(2*element[i].nod[j], 2*element[i].nod[l] + 1, 0);
                        Triplet <float> trpl21(2*element[i].nod[j] + 1, 2*element[i].nod[l], 0);
                        Triplet <float> trpl22(2*element[i].nod[j] + 1, 2*element[i].nod[l] + 1, 0);

                        tripl.push_back(trpl11);
                        tripl.push_back(trpl12);
                        tripl.push_back(trpl21);
                        tripl.push_back(trpl22);   
		}

		if ((i >= NE) and (i < NE + xBorder.size() / 2)) {	
			vector <int> infnod {NP + 2*(i - NE), xBorder[2*(i - NE)], xBorder[2*(i - NE) + 1], NP + 2*(i - NE) + 1};
			for (int j = 0; j < 4; j++)
				for (int l = 0; l < 4; l++) {
	
					Triplet <float> trpl11(2*infnod[j], 2*infnod[l], infel[i - NE].k(2*j, 2*l));
                	        	Triplet <float> trpl12(2*infnod[j], 2*infnod[l] + 1, infel[i - NE].k(2*j, 2*l + 1));
	                	        Triplet <float> trpl21(2*infnod[j] + 1, infnod[l], infel[i - NE].k(2*j + 1, 2*l));
	                        	Triplet <float> trpl22(2*infnod[j] + 1, 2*infnod[l] + 1, infel[i - NE].k(2*j + 1, 2*l + 1));

	        	                tripl.push_back(trpl11);
        	        	        tripl.push_back(trpl12);
                	        	tripl.push_back(trpl21);             	
					tripl.push_back(trpl22);
				}
		}
		
		if ((i >= NE + xBorder.size() / 2) and (i < NI)) {
	                vector <int> infnod {NP + 2*(i - NE), yBorder[2*(i - NE)], yBorder[2*(i - NE) + 1], NP + 2*(i - NE) + 1};
        	        for (int j = 0; j < 4; j++)
                	        for (int l = 0; l < 4; l++) {
	
        	                        Triplet <float> trpl11(2*infnod[j], 2*infnod[l], infel[i - NE].k(2*j, 2*l));
                	                Triplet <float> trpl12(2*infnod[j], 2*infnod[l] + 1, infel[i - NE].k(2*j, 2*l + 1));
                        	        Triplet <float> trpl21(2*infnod[j] + 1, infnod[l], infel[i - NE].k(2*j + 1, 2*l));
                                	Triplet <float> trpl22(2*infnod[j] + 1, 2*infnod[l] + 1, infel[i - NE].k(2*j + 1, 2*l + 1));

                                	tripl.push_back(trpl11);
	                                tripl.push_back(trpl12);
        	                        tripl.push_back(trpl21);
                	                tripl.push_back(trpl22);
                        	}
		}
        }


	K.setFromTriplets(tripl.begin(), tripl.end());
        K.makeCompressed();	

	VectorXf R (2*(NP + NI));
	for (int i = 0; i < 2*(NI + NP); i++)
		R(i) = 0;
	
	for (int i = 0; i < xBorder.size() / 2; i++) {
		if (nodesC[xBorder[2*i]][0] == minInColumn(nodesC, NP, 0)) {
			R(2*NP + 2*i) += -1 * F;
			R(2*NP + 2*i + 2) += -1 * F;
		}
		else {
			R(2*NP + 2*i) += F;
                        R(2*NP + 2*i + 2) += F;
		}
	}

	SimplicialLDLT<SparseMatrix <float> > solver(K);
	VectorXf U = solver.solve(R);
	vector <float> vU;
	for (int i = 0; i < U.size(); i++)
		vU.push_back(U(i));
	
	for (int i = 0; i < NE; i++) {
		for (int j = 0; j < 3; j++) {
			element[i].u(2*j) = vU[2*element[i].nod[j]];
			element[i].u(2*j + 1) = vU[2*element[i].nod[j] + 1];
		}
		element[i].sigma = D * element[i].B * element[i].u;
		element[i].epsilon = element[i].B * element[i].u;
		cout << element[i].epsilon << "\n------------\n";
	}
	
	for (int i = 0; i < xBorder.size() / 2; i++) {
		vector <int> infnod {NP + 2*i, xBorder[2*i], xBorder[2*i + 1], NP + 2*i + 1};
                for (int j = 0; j < 4; j++) {
                        infel[i].u(2*j) = vU[infnod[j]];
                        infel[i].u(2*j + 1) = vU[infnod[j] + 1];
                }
                infel[i].sigma = D * infel[i].B * infel[i].u;
		infel[i].epsilon = infel[i].B * infel[i].u;
		cout << infel[i].epsilon << "\n------------\n";

        }

	for (int i = xBorder.size() / 2; i < NI; i++) {
                vector <int> infnod {NP + 2*i, yBorder[2*(i - xBorder.size() / 2)], yBorder[2*(i - xBorder.size() / 2) + 1], NP + 2*i + 1};
                for (int j = 0; j < 4; j++) {
                        infel[i].u(2*j) = vU[infnod[j]];
                        infel[i].u(2*j + 1) = vU[infnod[j] + 1];
                }
		infel[i].sigma = D * infel[i].B * infel[i].u;
		infel[i].epsilon = infel[i].B * infel[i].u;
		cout << infel[i].epsilon << "\n------------\n";
        }


	return 0; 
}


