#include <cmath>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>
#include "gauss_kronrod.h"

static constexpr double MAX_DEPTH = 20;

const std::vector<double> nodes_gk = {
    0.000000000000000000000000000000000,
    0.148874338981631210884826001129720,
    0.294392862701460198131126603103866,
    0.433395394129247190799265943165784,
    0.562757134668604683339000099272694,
    0.679409568299024406234327365114874,
    0.780817726586416897063717578345042,
    0.865063366688984510732096688423493,
    0.930157491355708226001207180059508,
    0.973906528517171720077964012084452,
    0.995657163025808080735527280689003
};
const std::vector<double> weights_g = {
    0.295524224714752870173892994651338,
    0.269266719309996355091226921569469,
    0.219086362515982043995534934228163,
    0.149451349150580593145776339657697,
    0.066671344308688137593568809893332
};
const std::vector<double> weights_k = {
    0.149445554002916905664936468389821,
    0.147739104901338491374841515972068,
    0.142775938577060080797094273138717,
    0.134709217311473325928054001771707,
    0.123491976262065851077958109831074,
    0.109387158802297641899210590325805,
    0.093125454583697605535065465083366,
    0.075039674810919952767043140916190,
    0.054755896574351996031381300244580,
    0.032558162307964727478818972459390,
    0.011694638867371874278064396062192
};

/*funzione di singola integrazione. Restituisce integrale su intervallo specificato ed errore associato*/
std::pair<double, double> single_integral(double a, double b, const std::function<double(double)> &fnc) {
    double half_lenght = (b-a)/2;
    double center = (a+b)/2;
    double f_c = fnc(center);

    /*valutazione nodo origine nel kronrod (dispari) e inizializzazione a 0. nel gauss (pari)*/
    double res_gauss = 0.;
    double res_kron  = weights_k[0] * f_c;

    /*somma su nodi comuni, i pari del vector nodes_gk*/
    for (int i=1; i<nodes_gk.size(); i+=2){
        double temp = fnc(half_lenght*nodes_gk[i]+center) + fnc(-half_lenght*nodes_gk[i]+center);
        res_gauss += temp*weights_g[i/2];
        res_kron  += temp*weights_k[i];
    }
    /*somma su nodi mancanti, i dispari di nodes_gk*/
    for (int i=2; i<nodes_gk.size(); i+=2){
        double temp = fnc(half_lenght*nodes_gk[i]+center) + fnc(-half_lenght*nodes_gk[i]+center);
        res_kron += temp*weights_k[i];
    }

    res_gauss *= half_lenght;
    res_kron *= half_lenght;

    /*valutazione errore. Dettagli nel testo di tesi*/
    double err = pow(200*std::abs(res_gauss-res_kron), 1.5);

    return {res_kron, err};
}

/*funzione di integrazione ricorsiva che sfrutta la prima e ritorna valore dell'integrale con precisione richiesta*/
std::pair<double, double> integrate(double a, double b, const std::function<double(double)> &fnc, double eps_abs, int depth) {
    auto [res, err] = single_integral(a, b, fnc);
    if (err < eps_abs || depth > MAX_DEPTH) {
        return {res, err};
    }

    depth++;

    double midpoint = (a + b) * 0.5;
    auto [res_1, err_1] = integrate(a, midpoint, fnc, eps_abs, depth);
    auto [res_2, err_2] = integrate(midpoint, b, fnc, eps_abs, depth);

    return {res_1 + res_2, err_1 + err_2};
}
