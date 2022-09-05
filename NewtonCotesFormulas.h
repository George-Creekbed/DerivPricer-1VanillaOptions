#ifndef NEWTONCOTESFORMULAS_H_INCLUDED
#define NEWTONCOTESFORMULAS_H_INCLUDED

#include <iostream>
#include <functional>

namespace NewtonCotes {
    enum class Formula {
        Midpoint = 1,
        Trapezoidal,
        Simpsons
    };

    std::ostream& operator<< (std::ostream&, const Formula&);

    template<typename Function, typename Quadrature>
    double convergeToTol(const double, const double, int, const double, 
                         const Function&, const Quadrature&);

    template<typename Function>
    double printConvergenceValues(const Formula&, const double, const double, 
                                  const int, const double, const Function&);

    template<typename Function> double Midpoint(const double, const double, 
                                                const int, const Function&);

    template<typename Function> double Trapezoidal(const double, const double, 
                                                   const int, const Function&);

    template<typename Function> double Simpsons(const double, const double, 
                                                const int, const Function&);

    std::function<double (double, double, int, std::function<double (double) >) >
    createQuadrature(const Formula&);
}

#include "NewtonCotesFormulas.cpp"
#endif // NEWTONCOTESFORMULAS_H_INCLUDED
