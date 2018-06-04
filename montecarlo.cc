#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <montecarlo.h>
#include <gsl/gsl_randist.h>

using namespace std;

extern gsl_rng *rng;
#define sqr(x) ((x)*(x))

extern vector<neighbors> allnb;

//-----------------------------------------------------------------------------
double sigma(vector<double> &x)
{
	double result, sum = 0, a = avg(x);
	for (int i = 0; i < x.size(); i++)
	{
		sum += sqr(x[i]-a);
	}
	result = sqrt(sum/(x.size()-1));

	return result;
}
//-----------------------------------------------------------------------------
double avg(vector<double> &x)
{
	double s=0;
	for (int i = 0; i < x.size(); i++)
	{
		s += x[i];
	}
	return s/x.size();
}
//-----------------------------------------------------------------------------
void nb(vector<neighbors> &allnb, int L)  //Determine nearest neighbors of the site x (up,down, left, right)
{
	int V = L*L;
	for (int i=1; i<=V; i++)
	{	//up
		neighbors neigh(4);

		if ( (i-1) % L != 0)
		{
			neigh[0] = (i - 1);
		}
		else
		{
			neigh[0] = (i+L-1);
		}
		//down
		if (i % L == 0)
		{
			neigh[1] = (i-L+1);
		}
		else
		{
			neigh[1] = (i+1);
		}
		//left
		if (i <= L)
		{
			neigh[2] = (i+V-L);
		}
		else
		{
			neigh[2] = (i-L);
		}
		//right
		if (i > V - L)
		{
			neigh[3] = (i-V+L);
		}
		else
		{
			neigh[3] = (i+L);
		}
		allnb.push_back(neigh);
	}
	return;
}
//-----------------------------------------------------------------------------
void spin_config_hot(int L, vector<int> &spins)		//initialize a spin config
{

	int V = L*L;
	long seed2 = time(NULL);
	gsl_rng_set(rng,seed2);
	for (int i=0; i<V; i++)
		{
		if (gsl_rng_uniform(rng) < (1./2.))
		{
			spins.push_back(1);
		}
		else
		{
			spins.push_back(-1);
		}
	}
	return;
}

//-----------------------------------------------------------------------------
void spin_config_cold(int L, vector<int> &spins)		//initialize a spin config
{

	int V = L*L;
	for (int i=0; i<V; i++)
	{

		spins.push_back(1);
	}
	return;
}

//-----------------------------------------------------------------------------

int compE(int L, vector<int> &spins)  // Compute the energy of a config
{
	int E = 0, V = L*L;
	for (int i=0; i<V; i++)
	{
		int Esite = 0;
		for (int j=0; j<4; j+=2)
		{
			Esite = Esite - spins[i]*spins[allnb[i][j]]; 		// i-th site, j-th neighbor
		}
		E += Esite;
	}
	return E;
}

//-----------------------------------------------------------------------------

double compM(int L, vector<int> &spins)  //Compute magnetization of a config
{
	int V = L*L;
	double M = 0.;
	for (int i=0; i<V; i++)
	{
		M += spins[i];
	}
	return abs(M);
}
//-----------------------------------------------------------------------------

double comp_e(int L, vector<int> &spins) //Compute energy density of current config
{
	int V = L*L;
	double E = compE(L, spins);
	return E/V;
}

//-----------------------------------------------------------------------------

double comp_m(int L, vector<int> &spins)  //Compute magnetization density of a config
{
	int V = L*L;
	double M = abs(compM(L, spins));
	return M/V;
}

//-----------------------------------------------------------------------------

double lookup_exp(int q) 	//exp(beta*2*J*q), for beta = 1, J = 1
	{
		if (q == -4)
		{
			return 0.000335462627902511838821389125780861019310900133720319360;
		}
		if (q == -2)
		{
			return 0.018315638888734180293718021273241242211912067553475594769;
		}
		if (q == 0)
		{
			return 1.;
		}
		if (q == 2)
		{
			return 54.59815003314423907811026120286087840279073703861406872582;
		}
		else
		{
			return 2980.957987041728274743592099452888673755967939132835702208;
		}
	}

//-----------------------------------------------------------------------------

void flip(int L, vector<int> &spins, double T)
{
	//Perform a random single-spin flip on an Ising configs
	int V=L*L, flipped, nb_sum = 0;  // Site that will be flipped (j)
	double rho;
	flipped = gsl_rng_uniform_int(rng, V);

	for (int i=0; i<=3; i++)
	{
		nb_sum = nb_sum + spins[allnb[flipped][i]];  // Locall energy, nearest neighbors
	}

	int dE = -2*(-spins[flipped])*nb_sum; //Calculate the energy diff between the new and old configs

	if (dE <= 0)		// Metropolis accept-reject step
	{
		spins[flipped] = -spins[flipped];
	}
	else
	{
		//rho = lookup_exp(nb_sum);
		rho = exp(-dE/T);
		//rho = lookup_exp(-spins[flipped]*nb_sum);

		//cout << "rho= " << rho << " " << lookup_exp(-spins[flipped]*nb_sum) << endl;

		if (gsl_rng_uniform(rng) <= rho)
		{
			spins[flipped] = -spins[flipped];
		}
		else
		{
			return;
		}
	}
	return;
}

//-----------------------------------------------------------------------------

double spec_heat_dens(vector<double> &energies, int i_equil, double T, int L)
{// Calculates the specific heat density c=C/V, after equilibration
	//energies - vector containing E-s
	//i_equil - the index (MC time) of equilibration
	int V = L*L;
	vector<double> energies_equil;
	for (int j=i_equil; j<energies.size(); j++)
	{
		energies_equil.push_back(energies[j]);
	}
	double C = sqr(1./T)*sqr(sigma(energies_equil));
	return C/V;
}

//-----------------------------------------------------------------------------

double magn_suscept_dens(vector<double> &magnetiz, int i_equil, double T, int L)
{// Calculates the magnetic susceptibility density chi = chi_V/V, after equilibration
	//magnetiz - vector containing M-s
	//i_equil - the index (MC time) of equilibration
	int V = L*L;
	vector<double> magnetiz_equil;
	for (int j=i_equil; j<magnetiz.size(); j++)
	{
		magnetiz_equil.push_back(magnetiz[j]);
	}
	double chi = (1./T) * sqr(sigma(magnetiz_equil));
	return chi/V;
}

//-----------------------------------------------------------------------------

double autocovariance(vector<double> &y, int i_equil, int t)
{// y - primary quantity sample; i_equal - equilibration; t - time-shift
	// get the equilibrated vector:
	vector<double> y_equil;
	for (int j=i_equil; j<y.size(); j++)
	{
		y_equil.push_back(y[j]);
	}

	int N = y_equil.size();
	double pre = 1./(N-t);

	//first sum
	double first = 0.;
	double second = 0.;
	double third = 0.;

	for (int i=0; i<(N-t); i++)
	{
		first += y_equil[i+t];
		second += y_equil[i];
		third += y_equil[i]*y_equil[i+t];
	}

	return pre*third - (pre*second * pre*first);
}

//-----------------------------------------------------------------------------

double autocorrelation(vector<double> &y, int i_equil, int t)
{// y - primary quantity sample; i_equal - equilibration; t - time-shift
	double C_t = autocovariance(y, i_equil, t), C_0 = autocovariance(y, i_equil, 0);
	return C_t/C_0;
}

//-----------------------------------------------------------------------------

double int_auto_time(vector<double> &y, int i_equil) //Integrated autocorrelation time
{// y - primary quantity sample; i_equal - equilibration; t - time-shift
	double sum_rho = 0.5;
	for (int t=1; t<10000; t++)
	{
		double rho = autocorrelation(y, i_equil, t);
		if (rho > 0.)
		{
			sum_rho += rho;
		}
		else
		{
			break;
		}
	}
	return sum_rho;
}

//-----------------------------------------------------------------------------

double std_error_corr(vector<double> &y, int i_equil) // Std error on estimates of primary quantities <y>,  given by sample means y_bar of CORRELATED MC measurements
{// y - primary quantity sample; i_equal - equilibration; t - time-shift
	int N = y.size();
	return sqrt(2*int_auto_time(y, i_equil)/N) * sigma(y);
}

//-----------------------------------------------------------------------------

void effect_obs(vector <double> &g, vector<double> &x, vector<double> &y, double T, int L, int i_equil) //effective observable g(x,y) = f_x_bar * x + f_y_bar * y
{	// vector x := avg(E^2) ; y := avg(E)
	// Parameters: T - temperature, L - size, i_equil - equilibration time
	int V = L*L;
	double A = sqr(1./T)/V;
	vector<double> y_equil, x_equil;

	for (int j=i_equil; j<y.size(); j++)
	{
		y_equil.push_back(y[j]);
	}

	double avgy = avg(y_equil);

	for (int j=i_equil; j<x.size(); j++)
	{
		x_equil.push_back(x[j]);
	}

	for (int i = 1; i<x_equil.size(); i++)
	{
		g.push_back(A*x_equil[i] - 2*A*avgy*y_equil[i]);
	}

	return ;
}

//-----------------------------------------------------------------------------

double blocking(vector<double> &sample, int n_B, int i_equil, double T, int L, double (*f)(vector<double> &, int , double, int))
{
	vector<double> measure, sample_eq, Q;
	double Q_bar_bl, sum = 0;
	for (int j=i_equil; j<sample.size(); j++)
	{
		sample_eq.push_back(sample[j]);
	}

	int Ni = sample_eq.size()/n_B; //size of the block

	for (int j=0; j<n_B; j++)
	{
		Q.clear();
		for (int i=0; i<Ni; i++)
		{
			Q.push_back(sample_eq[j*Ni+i]);
		}
		//vector<double> &magnetiz, int i_equil, double T, int L
		measure.push_back((*f)(Q, i_equil, T, L));  // Q_(i) for the average (Q^bar_block)
	}
	Q_bar_bl = avg(measure);
	//Sum of the squares for std dev:
	for (int j=0; j<n_B; j++)
	{
		sum += sqr(measure[j] - Q_bar_bl);
	}

	return sqrt(sum/n_B)/sqrt(n_B);
}

//-----------------------------------------------------------------------------

double bootstrap(vector<double> &sample, int M, int i_equil, double T, int L, double (*f)(vector<double> &, int , double, int))
{// M - # of pseudo-samples
	vector<double> pseudo, measure, sample_eq;
	for (int j=i_equil; j<sample.size(); j++)
	{
		sample_eq.push_back(sample[j]);
	}

	int V = sample_eq.size();
	double Q_Boot, sum = 0;
	int autot = V/(2*int_auto_time(sample_eq, i_equil)); // Taking into account the autocorrelation time
	//cout << autot << endl;
	for (int j=0; j<M; j++) // Generate M - pseudo-samples
	{
		pseudo.clear();
		for (int i=0; i<autot; i++) // choose N/2*tau - configurations (because of autocorrelation)
		{
			pseudo.push_back(sample[gsl_rng_uniform_int(rng, V)]);
		}
		measure.push_back((*f)(pseudo, i_equil, T, L));  // Q_(i) for the average (Q^bar_Bootstrap)
	}

	Q_Boot = avg(measure);

	for (int j=0; j<M; j++)
	{
		sum += sqr(measure[j] - Q_Boot);
	}

	return sqrt(sum/M);
}
