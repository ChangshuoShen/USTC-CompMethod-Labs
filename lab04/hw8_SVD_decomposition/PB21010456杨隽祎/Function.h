#pragma once
#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
using namespace std;
void print(vector<vector<double>> A);
vector<vector<double>>transpose(vector<vector<double>> A);
double vec_infinite_norm(vector<double> z);
void house(vector<double>& v, double& beita, vector<double> x);
vector<vector<double>>  put(vector<vector<double>> A, int a, int b, int c, int d);
vector<double>  get(vector<vector<double>> A, int a, int b, int c, int d);
void copy(vector<vector<double> >& temp, vector<vector<double> >& A, int r, int c);
vector<vector<double>> unit(int n);
vector<double> mutilax(double a, vector<double>b);
vector<double> mutilAx(vector<vector<double>> A, vector<double> x);
vector<double> mutilxA(vector<vector<double>> A, vector<double> x);
vector<vector<double>> substract(vector<vector<double>> A, vector<vector<double>> B);
vector<vector<double> > mutil_vec(double beta, vector<double> u, vector<double> v);
void double_diagonal(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V);
int sign(double x);
void wilkinson(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V);
double mat_infinite_norm(vector<vector<double>> A);
vector<vector<double>> mutil(vector<vector<double>> m1, vector<vector<double>> m2);
int SVD(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V, vector<vector<double>>& C);
int Partition(vector<double>& A, int p, int r);
void QuickSort(vector<double>& A, int p, int r);
double  zuida(vector<vector<double>> A);
