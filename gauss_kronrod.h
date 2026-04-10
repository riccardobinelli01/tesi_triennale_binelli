#include <cmath>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

/*definizione funzione integrate. Se non specificato altrimenti, il parametro di profondità depth e di accuracy eps_abs vengono settati di default*/
std::pair<double, double> integrate(double a, double b, const std::function<double(double)> &fnc, double eps_abs = 1.0e-10, int depth = 0);
