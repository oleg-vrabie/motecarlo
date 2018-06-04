/*
 * statistics.h
 *
 *  Created on: 11.04.2018
 *      Author: vro29448
 */
#include <vector>
using namespace std;

typedef vector<int> neighbors;
double avg(vector<double> &x);
double sigma(vector<double> &x);
void nb(vector<neighbors> &all, int L);
void spin_config_hot(int L, vector<int> &spins);
void spin_config_cold(int L, vector<int> &spins);
int compE(int L, vector<int> &spins);
double compM(int L, vector<int> &spins);
double comp_e(int L, vector<int> &spins);
double comp_m(int L, vector<int> &spins);
void flip( int L, vector<int> &spins, double T);
double lookup_exp(int q);
double magn_suscept_dens(vector<double> &magnetiz, int i_equil, double T, int L);
double spec_heat_dens(vector<double> &energies, int i_equil, double T, int L);
double autocovariance(vector<double> &y, int i_equil, int t);
double autocorrelation(vector<double> &y, int i_equil, int t);
double int_auto_time(vector<double> &y, int i_equil);
double std_error_corr(vector<double> &y, int i_equil);
void effect_obs(vector <double> &g, vector<double> &x, vector<double> &y, double T, int L, int i_equil);
double blocking(vector<double> &sample, int n_B, int i_equil, double T, int L, double (*f)(vector<double> &, int , double, int));
double bootstrap(vector<double> &sample, int M, int i_equil, double T, int L, double (*f)(vector<double> &, int , double, int));
