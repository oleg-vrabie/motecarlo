#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <montecarlo.h>

#define sqr(x) ((x)*(x))

using namespace std;
gsl_rng *rng;

vector<neighbors> allnb;

int main()
{
	//"Random" part
	// gsl function to initialize the rng
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	// set seed based on current time
	long seed = time(NULL);
	gsl_rng_set(rng,seed);

	vector<int> spins;
	vector<double> upE, upM, upEsq, upMsq, upe, upm, geE, geM;

	int L=64;
	double T;
	nb(allnb, L);		// Determine nearest neighbors


// Ex 3:a) autocorrelation function rho(t)
/*
	T = 1.0;
	upe.clear();
	upm.clear();
	spin_config_cold(L, spins);
	// Metropolis sweep
	for (int j=0; j<10000; j++)
	{
		for (int i=0; i<(L*L); i++)		//One sweep
		{
			flip(L, spins, T);
		}
		cout << j << endl;
		upe.push_back(comp_e(L, spins));
		upm.push_back(comp_m(L, spins));
	}

	FILE * temp1c;
	temp1c = fopen("temp1c.txt","w");
	fprintf(temp1c, "#displacement t\trho_e(t)\trho_m(t) T=1.0, cold start\n");
	for (int t=1; t<10000; t++)
	{
		double rhoe = autocorrelation(upe, 200, t);
		double rhom = autocorrelation(upm, 200, t);
		cout << t << endl;
		if (rhom && rhoe > 0.)
		{
			fprintf(temp1c, "%d\t%f\t%f\n", t, rhoe, rhom);
		}
		else
		{
			break;
		}
	}
	fclose(temp1c);
	cout << "Integrated autocorrelation time for e: " << int_auto_time(upe, 200) << endl;
	cout << "Integrated autocorrelation time for m: " << int_auto_time(upm, 200) << endl;
	/*
	 *
	 *
	//-----------------------------------------------------------------
	T = 2.0;
	upe.clear();
	upm.clear();
	spin_config_cold(L, spins);
	// Metropolis sweep
	for (int j=0; j<10000; j++)
	{
		for (int i=0; i<(L*L); i++)		//One sweep
		{
			flip(L, spins, T);
		}
		cout << j << endl;
		upe.push_back(comp_e(L, spins));
		upm.push_back(comp_m(L, spins));
	}

	FILE * temp2c;
	temp2c = fopen("temp2c.txt","w");
	fprintf(temp2c, "#displacement t\trho_e(t)\trho_m(t) T=1.0, cold start\n");
	for (int t=1; t<10000; t++)
	{
		double rhoe = autocorrelation(upe, 200, t);
		double rhom = autocorrelation(upm, 200, t);
		cout << t << endl;
		if (rhom && rhoe > 0.)
		{
			fprintf(temp2c, "%d\t%f\t%f\n", t, rhoe, rhom);
		}
		else
		{
			break;
		}
	}
	fclose(temp2c);
	//cout <<  spec_heat_dens(upE, 100, T, L) << endl;
	//cout <<  magn_suscept_dens(upM, 100, T, L) << endl;
	//cout << compE(L, spins)<<' ' << compM(L, spins)<< '\n';
	cout << "Integrated autocorrelation time for e: " << int_auto_time(upe, 200) << endl;
	cout << "Integrated autocorrelation time for m: " << int_auto_time(upm, 200) << endl;
	*/

	//-----------------------------------------------------------------
	for (int T=1; T<4; T++)
	{
		geE.clear();
		geM.clear();
		upE.clear();
		upM.clear();
		upEsq.clear();
		upMsq.clear();
		spin_config_cold(L, spins);
		// Metropolis sweep

		for (int j=0; j<10000; j++)
		{
			for (int i=0; i<(L*L); i++)		//One sweep
			{
				flip(L, spins, T);
			}
			//cout << j << endl;
			//upe.push_back(comp_e(L, spins));
			//upm.push_back(comp_m(L, spins));
			upE.push_back(compE(L, spins));
			upM.push_back(compM(L, spins));
			upEsq.push_back(sqr(compE(L, spins)));
			upMsq.push_back(sqr(compM(L, spins)));
		}
		/*
		FILE * temp3c;
		temp3c = fopen("temp3c.txt","w");
		fprintf(temp3c, "#displacement t\trho_e(t)\trho_m(t) T=1.0, cold start\n");
		for (int t=1; t<10000; t++)
		{
			double rhoe = autocorrelation(upE, 200, t);
			double rhom = autocorrelation(upM, 200, t);
			cout << t << endl;
			if (rhom && rhoe > 0.)
			{
				fprintf(temp3c, "%d\t%f\t%f\n", t, rhoe, rhom);
			}
			else
			{
				break;
			}
		}
		fclose(temp3c);

		*/
		effect_obs(geE, upEsq, upE, T, L, 200);
		effect_obs(geM, upMsq, upE, T, L, 200);
		//cout <<  spec_heat_dens(upE, 100, T, L) << endl;
		//cout <<  magn_suscept_dens(upM, 100, T, L) << endl;
		//cout << compE(L, spins)<<' ' << compM(L, spins)<< '\n';
		//cout << "Integrated autocorrelation time for e: " << int_auto_time(upe, 200) << endl;
		cout << "Temperature T=" <<T << endl;
		cout << "Specific heat density, c: " << spec_heat_dens(upE, 200, T, L) << endl;
		cout << "Magnetic susceptibility density, chi: " << magn_suscept_dens(upM, 200, T, L) << endl;
		//cout << "Integrated autocorrelation time for m: " << int_auto_time(upm, 200) << endl;
		cout << "Std error of c using error propagation: " << std_error_corr(geE, 200) << endl;
		cout << "Std error of c using blocking method: " << blocking(upE, 20, 200, T, L, spec_heat_dens) << endl;
		cout << "Std error of c using bootstrap method: " << bootstrap(upE, 1000, 200, T, L, spec_heat_dens) << endl;

		cout << "Std error of chi using error propagation : " << std_error_corr(geM, 200) << endl;
		cout << "Std error of chi using blocking method: " << blocking(upM, 20, 200, T, L, magn_suscept_dens) << endl;
		cout << "Std error of chi using bootstrap method: " << bootstrap(upM, 1000, 200, T, L, magn_suscept_dens) << endl;
		cout << endl;
	}

	gsl_rng_free(rng);
	return 0;
}
