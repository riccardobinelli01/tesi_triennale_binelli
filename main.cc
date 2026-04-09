#include <cmath>
#include <functional>
#include <iostream>
#include <utility> /*std::pair*/
#include <vector>
#include "gauss_kronrod.h"
#include <iomanip> /*std::setprecision(int numero_cifre_significative)*/

/*routine che sfrutta gauss_kronrod 10-21 per calcolare fattore di Gamow e il tempo di vita medio degli isotopi dell'uranio. Usa il sistema di unità naturali*/

const double charge_squared = 197.32/137.04; //e^2 scritto in MeV*fm
//const double atomic_number = 92;

/*potenziale restituisce in MeV, vanno passati sia r che w in fermi (fm)*/
double potential(double r, double w, double V_0, double Z) {
    double temp = 2*(Z - 2)*charge_squared; //2(Z-2)e^2 * (fattore di conversione 1/fm -> MeV)

    //diviso in due regioni: dentro e fuori dal nucleo
    if(r<=w) return temp*((3*w*w-r*r)/(2*w*w*w)) + V_0/(exp((r-w)/0.8)+1); //0.8 corrisponde ad "a" in fermi (fm)
    else return temp/r + V_0/(exp((r-w)/0.8)+1);
}

int main()
{
    double mass = 3727.38; //MeV
    double atomic_number = 92;
    std::vector<double> mass_numbers = {234,     232,      230,      228,     226,      224,      222,      220     };
    std::vector<double> widths =       {8.01,    7.99,     7.96,     7.94,    7.92,     7.89,     7.87,     7.85    }; //fm
    std::vector<double> V_0 =          {-115.72, -115.85,  -116.00,  -115.78, -115.53,  -114.98,  -114.25,  -113.54 }; //MeV
    std::vector<double> E_alpha =      {4.27,    4.57,     4.86,     5.41,    5.99,     6.80,     7.72,     8.62    }; //MeV
    std::vector<double> P =            {0.22,    0.15,     0.12,     0.11,    0.10,     0.071,    0.044,    0.027   };
    //trovati con desmos
    /*std::vector<double> r_1 =          {9.07039, 9.06247,  9.04346,  9.04527, 9.04883,  9.05042,  9.06928,  9.08925 }; //fm
    std::vector<double> r_2 =          {60.7149, 56.71267, 53.32858, 47.907,  43.26826, 38.11425, 33.57214, 30.06692}; //fm*/
    //valori Pasquini
    std::vector<double> r_1 =          {9.07039, 9.05978,  9.04883,  9.04699, 9.04673,  9.05597,  9.07066,  9.08624 }; //fm
    std::vector<double> r_2 =          {60.7149, 56.7292,  53.3441,  47.9210, 43.2809,  38.1254,  33.5819,  30.0757 }; //fm

    std::vector<double> tempi_dimezz_bnlgov = {1.40902848e17, 7.3857312e14, 7742088000000, 2172830400, 1797120, 546, 0.268, 0.000396};
    std::vector<double> tempi_dimezz_articolo = {1.41e17, 7.39e14, 7.73e12, 2.17e9, 1.80e6, 3.28e4, 3.5e-1, 9e-4, 1e-6};

    //stampa della tabellona con parametri
    /*
    for(int i=0; i<mass_numbers.size(); i++){
        std::cout << "U" << mass_numbers[i]+4 << " | " << mass_numbers[i] << "\t" << widths[i] << "\t" << V_0[i] << "\t" << E_alpha[i] << "\t" << P[i] << "\t" << r_1[i] << "\t" << r_2[i] << std::endl;
    }
    std::cout << std::endl;
    */

    /*stampa in formato ristretto*/
    std::cout << "isotopo\t|\t#massa\tw\tV0\tE_alpha\tP\tR\tR_1\t|\ttau_1/2" << std::endl;
    for(int i=0; i<91; i++) std::cout<<"-";
    std::cout << std::endl;
    /*std::cout << "isotopo\t|\tP" << std::endl;
    for(int i=0; i<26; i++) std::cout<<"-";
    std::cout << std::endl;*/

    for(int i=0; i<mass_numbers.size(); i++){
        //costruzione funzione
        double particle_energy = E_alpha[i];
        double width = widths[i];
        //double width = 1.3*std::pow(mass_numbers[i], 1./3.);
        double V = V_0[i];
        std::function<double(double)> fnc = [particle_energy, width, V, atomic_number](double x) -> double {
            return std::pow(std::abs(particle_energy - potential(x, width, V, atomic_number)), 0.5);
        };

        //integrazione tramite integrate
        double ingresso = r_1[i];
        double uscita = r_2[i];
        auto [gamow, err] = integrate(ingresso, uscita, fnc);
        gamow *= sqrt(2.*mass)/197.32; //ultimo termine è hbar in MeV*fm

        //fattore hit della barriera e preformazione P
        double epsilon = E_alpha[i]-(V_0[i] + 1.5*2.*(atomic_number - 2.)*charge_squared/widths[i]);
        //std::cout << epsilon << std::endl;
        double hit_factor = sqrt((2.*epsilon)/mass)*(1./(2.*ingresso));
        //std::cout << hit_factor << std::endl;

        double lambda_BP = hit_factor*exp(-2.*gamow);
        double lambda_finale = lambda_BP*P[i]*0.5;
            /*simulazione di P
            double lambda_correct = 0.5*lambda_BP;
            double lambda_exp = (0.6931471806*3.3356e-24)/tempi_dimezz_bnlgov[i]; //primo termine è ln(2) e conversione secondi -> fm
            double lambda_finale = lambda_BP*0.5*(lambda_exp/lambda_correct);*/

        /*stampa per il terminale
        std::cout << "Isotopo di Uranio U" << mass_numbers[i]+4 << std::endl;
        std::cout << "Fattore di Gamow:          " << gamow << std::endl;
        std::cout << "Prob decad lambda_BP:      " << lambda_BP << std::endl;
        std::cout << "Prob decad lambda_correct: " << 0.5*lambda_BP << std::endl;
        std::cout << "Prob decad lambda_finale:     " << lambda_finale << std::endl;
        std::cout << "Tempo di vita medio:       " << 3.3356e-24/lambda_finale << std::endl; //ultimo termine per convertire fm -> secondi
        std::cout << "Tempo di dimezzamento:     " << (0.6931471806*3.3356e-24)/lambda_finale << std::endl;
        std::cout << std::endl;*/

        //stampa in formato ristretto
        std::cout << "U" << mass_numbers[i]+4 << "\t|\t" << mass_numbers[i] << "\t" << width << "\t" << V_0[i] << "\t" << E_alpha[i] << "\t" << P[i] << "\t" << r_1[i] << "\t" << r_2[i] << "\t|\t" << (0.6931471806*3.3356e-24)/lambda_finale << std::endl;
        //std::cout << tempi_dimezz_bnlgov[i] << " & " << (0.6931471806*3.3356e-24)/lambda_finale << "\\" << "\\" << std::endl;

        /*simulazione di P
        double lambda_correct = 0.5*lambda_BP;
        double lambda_exp = (0.6931471806*3.3356e-24)/tempi_dimezz_articolo[i]; //primo termine è ln(2) e conversione secondi -> fm
        std::cout << "U" << mass_numbers[i]+4 << "\t|\t" << lambda_exp/lambda_correct << std::endl;*/

    }
    //for(int i=0; i<tempi_bnlgov.size(); i++) std::cout << "U" << mass_numbers[i]+4 << " | " << tempi_dimezz_articolo[i] << std::endl;

    //per stampare i vettori da mettere in Colab
    //std::cout << std::endl;

    return 0;
}
