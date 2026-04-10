#include <cmath>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>
#include "gauss_kronrod.h"
#include <iomanip>

/*Routine che sfrutta regola di gauss_kronrod 10-21 per calcolare il tempo di dimezzamento degli isotopi pari dell'Uranio, da 224U a 238U. Usa il sistema di unità naturali*/

const double charge_squared = 197.32/137.04; //e^2 scritto in MeV*fm

/*definizione potenziale elettromagnetico + nucleare. Restituisce in MeV, r e w sono in fm*/
double potential(double r, double w, double V_0, double Z) {
    double temp = 2*(Z - 2)*charge_squared;

    if(r<=w) return temp*((3*w*w-r*r)/(2*w*w*w)) + V_0/(exp((r-w)/0.8)+1);
    else return temp/r + V_0/(exp((r-w)/0.8)+1);
}

int main()
{
    double mass = 3727.38; /*massa particella alpha in MeV*/
    double atomic_number = 92;
    std::vector<double> mass_numbers = {234,     232,      230,      228,     226,      224,      222,      220     };
    std::vector<double> widths =       {8.01,    7.99,     7.96,     7.94,    7.92,     7.89,     7.87,     7.85    }; //fm
    std::vector<double> V_0 =          {-115.72, -115.85,  -116.00,  -115.78, -115.53,  -114.98,  -114.25,  -113.54 }; //MeV
    std::vector<double> E_alpha =      {4.27,    4.57,     4.86,     5.41,    5.99,     6.80,     7.72,     8.62    }; //MeV
    std::vector<double> P =            {0.22,    0.15,     0.12,     0.11,    0.10,     0.071,    0.044,    0.027   };
    std::vector<double> r_1 =          {9.07039, 9.05978,  9.04883,  9.04699, 9.04673,  9.05597,  9.07066,  9.08624 }; //fm
    std::vector<double> r_2 =          {60.7149, 56.7292,  53.3441,  47.9210, 43.2809,  38.1254,  33.5819,  30.0757 }; //fm

    std::cout << "isotopo\t|\t#massa\tw\tV0\tE_alpha\tP\tR\tR_1\t|\ttau_1/2" << std::endl;
    for(int i=0; i<91; i++) std::cout<<"-";
    std::cout << std::endl;

    for(int i=0; i<mass_numbers.size(); i++){
        /*costruzione funzione*/
        double particle_energy = E_alpha[i];
        double width = widths[i];
        double V = V_0[i];
        std::function<double(double)> fnc = [particle_energy, width, V, atomic_number](double x) -> double {
            return std::pow(std::abs(particle_energy - potential(x, width, V, atomic_number)), 0.5);
        };

        /*definizione intervallo e integrazione tramite integrate*/
        double ingresso = r_1[i];
        double uscita = r_2[i];
        auto [gamow, err] = integrate(ingresso, uscita, fnc);
        gamow *= sqrt(2.*mass)/197.32; //ultimo termine è hbar in MeV*fm

        /*calcolo costante di decadimento con fattore hit della barriera, esponenziale di Gamow e preformazione P*/
        double epsilon = E_alpha[i]-(V_0[i] + 1.5*2.*(atomic_number - 2.)*charge_squared/widths[i]);
        double hit_factor = sqrt((2.*epsilon)/mass)*(1./(2.*ingresso));
        double lambda_BP = hit_factor*exp(-2.*gamow);
        double lambda_finale = lambda_BP*P[i]*0.5;

        /*stampa tabella risultati. L'ultimo termine contiene la conversione fm -> s*/
        std::cout << "U" << mass_numbers[i]+4 << "\t|\t" << mass_numbers[i] << "\t" << width << "\t" << V_0[i] << "\t" << E_alpha[i] << "\t" << P[i] << "\t" << r_1[i] << "\t" << r_2[i] << "\t|\t" << (0.6931471806*3.3356e-24)/lambda_finale << std::endl;

    }

    return 0;
}
