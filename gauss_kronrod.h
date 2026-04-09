#include <cmath>
#include <functional>
#include <iostream>
#include <utility> /*std::pair*/
#include <vector>

/*se non specifico l'eps_abs e la depth (di partenza) li setto di default*/
std::pair<double, double> integrate(double a, double b, const std::function<double(double)> &fnc, double eps_abs = 1.0e-10, int depth = 0);
