#include"Function.h"
#include <cmath>
#include <iomanip> 
#include<vector>
#include<iostream>
#include<cstdlib>
#include<Windows.h>
#include<ctime>
#include<complex>
using namespace std;
void print(vector<vector<double>> A) 
{
	int m = A.size();
	int n = A[0].size();
	int i, j;
	for (i = 0; i < m; i++) 
	{
		cout << endl;
		for (j = 0; j < n; j++) 
		{
			cout << A[i][j] << "       ";
		}
	}

	cout << endl;
	cout << endl;
}

vector<vector<double>>transpose(vector<vector<double>> A) //将A转置
{
	
	int n = A.size();
	int i, j;
	double temp;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
	return A;
}
double vec_infinite_norm(vector<double> z)//向量的∞范数
{
	double infinite_norm = 0;

	for (int i = 0; i < z.size(); i++)
		if (fabs(z[i]) > infinite_norm)
		{
			infinite_norm = fabs(z[i]);

		}
	return infinite_norm;
}
void house(vector<double>& v, double& beita, vector<double> x)//构造house
{
	double a;
	int n = x.size();
	double nita = vec_infinite_norm(x);
	for (int i = 0; i < n; i++) {
		x[i] = x[i] / nita;
	}
	double sigema = 0;
	for (int i = 1; i < n; i++) {
		sigema = sigema + x[i] * x[i];
	}

	for (int i = 1; i < n; i++) {
		v[i] = x[i];
	}
	if (sigema == 0)
		beita = 0;
	else {
		a = sqrt(x[0] * x[0] + sigema);
		if (x[0] <= 0) {
			v[0] = x[0] - a;
		}
		else {
			v[0] = -sigema / (x[0] + a);
		}
		beita = 2 * v[0] * v[0] / (sigema + v[0] * v[0]);
		double temp = v[0];
		for (int i = 0; i < n; i++) {
			v[i] = v[i] / temp;
		}
	}

}
vector<vector<double>>  put(vector<vector<double>> A, int a, int b, int c, int d) //带入
{
	int i, j;
	int m = b - a + 1;
	int n = d - c + 1;
	vector<vector<double>> x(m, vector<double>(n));
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			x[i][j] = A[a + i][c + j];
		}
	}
	return x;
}
vector<double>  get(vector<vector<double>> A, int a, int b, int c, int d) {
	//a- b 行 c -d 列
	//取出 
	int i;
	if (c == d) {
		//只取出一列
		int m = b - a + 1;
		vector<double> x(m);
		for (i = 0; i < m; i++) {
			x[i] = A[a + i][c];
		}
		return x;
	}
	else if (a == b) {
		//一行
		int m = d - c + 1;
		vector<double> x(m);
		for (i = 0; i < m; i++) {
			x[i] = A[a][c + i];
		}
		return x;
	}
	else {
		cout << " error get" << endl;
	}
}
void copy(vector<vector<double> >& temp, vector<vector<double> >& A, int r, int c)//复制
{
	int m = temp.size();
	int n = temp[0].size();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[r + i][c + j] = temp[i][j];
}
vector<vector<double>> unit(int n) 
{
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0; i < n; i++) 
	{
		A[i][i] = 1;
	}
	return A;
}
vector<double> mutilax(double a, vector<double>b) //a乘以向量
{
	int n = b.size();
	vector<double>c(n);
	for (int i = 0; i < n; i++) {
		c[i] = a * b[i];
	}
	return c;
}
vector<double> mutilAx(vector<vector<double>> A, vector<double> x)//矩阵向量相乘
{

	int m = A.size();
	int n = A[0].size();
	vector<double>y(m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			y[i] = y[i] + A[i][j] * x[j];
		}
	}
	return y;
}
vector<double> mutilxA(vector<vector<double>> A, vector<double> x)//向量矩阵相乘
{

	int m = A.size();
	int n = A[0].size();
	vector<double>y(n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			y[j] = y[j] + A[i][j] * x[i];
		}
	}
	return y;
}
vector<vector<double>> substract(vector<vector<double>> A, vector<vector<double>> B)
{
	int m = A.size();
	int n = A[0].size();
	vector<vector<double>> C(m, vector<double>(n));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
}
vector<vector<double> > mutil_vec(double beta, vector<double> u, vector<double> v) {
	int m = u.size();
	int n = v.size();
	vector<vector<double> > temp(m, vector<double>(n));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			temp[i][j] = u[i] * v[j] * beta;
	return temp;
}
void double_diagonal(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V) {
	int m = A.size();
	int n = A[0].size();
	vector<double>b(n);
	vector<double>c(n);
	vector<vector<double>> P(m, vector<double>(m));
	vector<vector<double>> Q(n, vector<double>(n));
	P = unit(m);
	Q = unit(n);
	for (int k = 0; k < n; k++) {
		double beta;
		vector<double>v(m - k);

		house(v, beta, get(A, k, m - 1, k, k));
		vector<double>u(n - k); u = mutilax(beta, mutilxA(put(A, k, m - 1, k, n - 1), v));
		vector<vector<double>> A1(m - k, vector<double>(n - k));
		A1 = substract(put(A, k, m - 1, k, n - 1), mutil_vec(1, v, u));
		copy(A1, A, k, k);
		for (int i = 0; i < m - 1 - k; i++) {
			A[k + 1 + i][k] = v[1 + i];
		}
		b[k] = beta;
		vector<double>p(m); p = mutilax(beta, mutilxA(put(P, k, m - 1, 0, m - 1), v));
		vector<vector<double>> P1(m - k, vector<double>(m));
		P1 = substract(put(P, k, m - 1, 0, m - 1), mutil_vec(1, v, p));
		copy(P1, P, k, 0);
		if (k < n - 2) {
			double beta1;
			vector<double>v1(n - 1 - k);
			house(v1, beta1, get(A, k, k, k + 1, n - 1));
			vector<double>u1(m - k); u1 = mutilax(beta1, mutilAx(put(A, k, m - 1, k + 1, n - 1), v1));
			vector<vector<double>> A2(m - k, vector<double>(n - k - 1));
			A2 = substract(put(A, k, m - 1, k + 1, n - 1), mutil_vec(1, u1, v1));
			copy(A2, A, k, k + 1);
			for (int i = 0; i < n - k - 2; i++) {
				A[k][k + 2 + i] = v1[1 + i];
			}
			c[k] = beta1;
			vector<double>q1(n); q1 = mutilax(beta1, mutilAx(put(Q, 0, n - 1, k + 1, n - 1), v1));
			vector<vector<double>> Q1(n, vector<double>(n - k - 1));
			Q1 = substract(put(Q, 0, n - 1, k + 1, n - 1), mutil_vec(1, q1, v1));
			copy(Q1, Q, 0, k + 1);
		}
	}
	U = P;
	V = Q;
}

int sign(double x) 
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

void wilkinson(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V) 
{
	int n = B.size();
	vector<vector<double>> P(n, vector<double>(n));
	P = unit(n);
	vector<vector<double>> Q(n, vector<double>(n));
	Q = unit(n);
	vector<double> delta(n + 1);
	vector<double> gamma(n);
	vector<double> Qc(n);
	vector<double> Qs(n);
	vector<double> Pc(n);
	vector<double> Ps(n);
	double temp3, temp4;
	for (int i = 1; i < n + 1; i++) 
	{
		delta[i] = B[i - 1][i - 1];
	}
	for (int i = 1; i < n; i++) 
	{

		gamma[i] = B[i - 1][i];
	}
	double afa = delta[n] * delta[n] + gamma[n - 1] * gamma[n - 1];
	double derta1 = (delta[n - 1] * delta[n - 1] + gamma[n - 2] * gamma[n - 2] - afa) / 2;
	double beta = delta[n - 1] * gamma[n - 1];
	double miu = afa - beta * beta / (derta1 + sign(derta1) * sqrt(afa * afa + beta * beta));
	double y = delta[1] * delta[1] - miu;
	double z = delta[1] * gamma[1];
	double c, s, temp1, temp2;
	for (int k = 1; k < n; k++) 
	{
		c = y / sqrt(y * y + z * z);
		s = -z / sqrt(y * y + z * z);
		gamma[k - 1] = sqrt(y * y + z * z);
		Qc[k] = c; Qs[k] = s;
		for (int i = 0; i < n; i++) 
		{
			temp3 = c * Q[i][k - 1] - s * Q[i][k];
			temp4 = s * Q[i][k - 1] + c * Q[i][k];
			Q[i][k - 1] = temp3;
			Q[i][k] = temp4;
		}
		y = c * delta[k] - s * gamma[k];
		z = -s * delta[k + 1];
		temp1 = s * delta[k] + c * gamma[k];
		temp2 = c * delta[k + 1];
		gamma[k] = temp1;
		delta[k + 1] = temp2;
		delta[k] = sqrt(y * y + z * z);
		c = y / sqrt(y * y + z * z);
		s = -z / sqrt(y * y + z * z);
		Pc[k] = c;
		Ps[k] = s;
		for (int i = 0; i < n; i++) {
			temp3 = c * P[k - 1][i] - s * P[k][i];
			temp4 = s * P[k - 1][i] + c * P[k][i];
			P[k - 1][i] = temp3;
			P[k][i] = temp4;
		}
		if (k < n - 1) {
			y = c * gamma[k] - s * delta[k + 1];
			z = -s * gamma[k + 1];
			temp1 = s * gamma[k] + c * delta[k + 1];
			temp2 = c * gamma[k + 1];
			delta[k + 1] = temp1;
			gamma[k + 1] = temp2;
		}
		else {
			temp1 = c * gamma[k] - s * delta[k + 1];
			temp2 = s * gamma[k] + c * delta[k + 1];
			gamma[k] = temp1;
			delta[k + 1] = temp2;
		}
	}
	for (int i = 1; i < n + 1; i++) 
	{
		B[i - 1][i - 1] = delta[i];

	}
	for (int i = 1; i < n; i++) 
	{

		B[i - 1][i] = gamma[i];
	}
	U = P;
	V = Q;
}

double mat_infinite_norm(vector<vector<double>> A)//矩阵无穷范数 
{

	int i, j, n;
	double max, sum;
	n = A.size();
	max = 0;
	for (i = 0; i < n; i++) {
		sum = 0;
		for (j = 0; j < n; j++) {
			sum = sum + abs(A[i][j]);
		}
		if (sum > max) {
			max = sum;
		}
	}

	return max;
}
vector<vector<double>> mutil(vector<vector<double>> m1, vector<vector<double>> m2)//两矩阵相乘
{

	int m = m1.size();
	int n = m1[0].size();
	int p = m2[0].size();
	vector<vector<double>> array;
	vector<double> temparay;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += m1[i][k] * m2[k][j];
			}
			temparay.push_back(sum);
		}
		array.push_back(temparay);
		temparay.erase(temparay.begin(), temparay.end());
	}
	return array;
}
int SVD(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V, vector<vector<double>>& C)
{

	int m = A.size();
	int n = A[0].size();
	vector<vector<double>> P(m, vector<double>(m));
	vector<vector<double>> Q(n, vector<double>(n));

	P = unit(m);
	Q = unit(n);

	double_diagonal(A, P, Q);
	double acc = 1e-7;
	vector<vector<double>> B(n, vector<double>(n));
	for (int i = 0; i < n - 1; i++) {

		B[i][i] = A[i][i];
		B[i][i + 1] = A[i][i + 1];
	}
	B[n - 1][n - 1] = A[n - 1][n - 1];
	int j = 0;
	while (1) {
		j++;
		for (int i = 0; i < n - 1; i++) {
			if (abs(B[i][i + 1]) < (abs(B[i][i]) + abs(B[i + 1][i + 1])) * acc)
				B[i][i + 1] = 0;
			if (abs(B[i][i]) <= mat_infinite_norm(B) * acc) { B[i][i] = 0; }
		}
		int q;
		for (q = n - 1; q > 0; ) {
			if (B[q - 1][q] == 0)q--;
			else { break; }
		}
		
		if (q == 0)break;
		int p;
		for (p = q; p > 0;) {
			if (B[p - 1][p] != 0) {
				p--; 
			}

			else { break; }
		}
		
		for (int i = p; i <= q; i++) {
			if (A[i][i] == 0) {
				cout << "对角线0";

				break;
			}
		}
		vector<vector<double>> inP(m, vector<double>(m));
		vector<vector<double>> inQ(n, vector<double>(n));
		inP = unit(m);
		inQ = unit(n);
		vector<vector<double>> B1(q - p + 1, vector<double>(q - p + 1));
		vector<vector<double>> Q1(q - p + 1, vector<double>(q - p + 1));
		vector<vector<double>> P1(q - p + 1, vector<double>(q - p + 1));
		B1 = put(B, p, q, p, q);
		wilkinson(B1, P1, Q1);
		copy(B1, B, p, p);
		copy(P1, inP, p, p);
		copy(Q1, inQ, p, p);
		P = mutil(inP, P);
		Q = mutil(Q, inQ);
	}
	U = P;
	V = Q;
	C = B;
	return j;
}

int Partition(vector<double>& A, int p, int r) {

	double x;
	int i, j;

	x = A[r];
	i = p - 1;

	for (j = p; j < r; j++) 
	{
		if (A[j] <= x) 
		{
			i++;
			swap(A[i], A[j]);
		}
	}
	swap(A[i + 1], A[r]);
	return (i + 1);
}

void QuickSort(vector<double>& A, int p, int r) //快排
{
	int q;
	if (p < r) 
	{
		q = Partition(A, p, r);
		QuickSort(A, p, q - 1);
		QuickSort(A, q + 1, r);
	}
}

double  zuida(vector<vector<double>> A) {
	int m = A.size();
	int n = A[0].size();
	int i, j;
	double temp = 0;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			if (temp < abs(A[i][j])) {
				temp = abs(A[i][j]);
			}
		}
	}

	return  temp;
}