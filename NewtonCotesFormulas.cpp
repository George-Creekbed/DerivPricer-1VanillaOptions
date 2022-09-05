#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include "NewtonCotesFormulas.h"

using namespace NewtonCotes;
using std::string;          using std::cout;
using std::setprecision;    using std::fixed;
using std::scientific;      using std::ios_base;
using std::abs;             using std::function;
using std::runtime_error;   using std::ostream;
using std::to_string;

// this prints the enum on cout
ostream& operator<< (ostream& os, const Formula& method) {
    if (method == Formula::Midpoint)
        os << "Midpoint";
    else if (method == Formula::Trapezoidal)
        os << "Trapezoidal";
    else if (method == Formula::Simpsons)
        os << "Simpsons";

    return os;
}

template<typename Function, typename Quadrature>
double convergeToTol(const double inf, const double sup, int intervals, const double tol,
                     const Function& func, const Quadrature& my_method) {
    double old_integral = my_method(inf, sup, intervals, func);
    intervals *= 2;
    double new_integral = my_method(inf, sup, intervals, func);
    double difference = abs(new_integral - old_integral);

    cout.setf(ios_base::scientific);
    cout.precision(8);
    while(difference > tol) {
        cout << string(9 - to_string(intervals).size(), ' ') << intervals <<
                string(3,' ') << fixed << new_integral << string(3,' ') <<
                scientific << difference << '\n';
        old_integral = new_integral;
        intervals *= 2;
        new_integral = my_method(inf, sup, intervals, func);
        difference = abs(new_integral - old_integral);
    }
    cout << fixed << string(9 - to_string(intervals).size(), ' ') <<
            intervals << string(3,' ') << new_integral << string(3,' ') <<
            scientific << difference << '\n';

    return new_integral;
}

template<typename Function>
double printConvergenceValues(const Formula& choice, const double inf, const double sup, const int intervals,
                              const double tol, const Function& func)
{
    function<double (double, double, int, function<double (double) >) >
        my_method = createQuadrature(choice);

    cout << choice << " method, with tolerance " << tol << "\n";
    cout << "intervals" << string(5,' ') << "integral " << string(9,' ') << "tol\n";
    cout << setprecision(8) << fixed << string(9 - to_string(intervals).size(), ' ') <<
            intervals << string(3,' ') << my_method(inf, sup, intervals, func) <<
            string(3,' ') << 0.0 << '\n';

    double result = convergeToTol(inf, sup, intervals, tol, func, my_method);

    return result;
}

template <typename Function>
double Midpoint(const double inf, const double sup, const int intervals, const Function& func) {
    double interval_width = (sup - inf) / intervals;
    double result = 0;
    for (int i = 1; i <= intervals; ++i)
        result += func(inf + (i - 0.5) * interval_width);

    return result * interval_width;
}

template <typename Function>
double Trapezoidal(const double inf, const double sup, const int intervals, const Function& func) {
    double interval_width = (sup - inf) / intervals;
    double result = (func(inf) + func(sup)) / 2.;
    for (int i = 1; i <= intervals - 1; ++i)
        result += func(inf + i * interval_width);

    return result * interval_width;
}

template <typename Function>
double Simpsons(const double inf, const double sup, const int intervals, const Function& func) {
    double interval_width = (sup - inf) / intervals;
    double result = (func(inf) + func(sup)) / 6;
    for (int i = 1; i <= intervals - 1; ++i)
        result += func(inf + i * interval_width) / 3 + 2 * func(inf + (i - 0.5) * interval_width) / 3;
    result += 2 * func(inf + (intervals - 0.5) * interval_width) / 3;

    return result * interval_width;
}

function<double (double, double, int, function<double (double) >) >
createQuadrature(const Formula& choice) {
    if (choice == Formula::Midpoint)
        return NewtonCotes::Midpoint<function<double (double)> >;
    else if (choice == Formula::Trapezoidal)
        return NewtonCotes::Trapezoidal<function<double (double)> >;
    else if (choice == Formula::Simpsons)
        return NewtonCotes::Simpsons<function<double (double)> >;
    else
        throw runtime_error("Pick between choice = 1, 2, or 3");

    return 0;  // never gets here
}
