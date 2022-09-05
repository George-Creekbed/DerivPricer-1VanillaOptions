#define _USE_MATH_DEFINES
/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG1E    - log2(e)
 * M_LOG9E    - log10(e)
 * M_LN1      - ln(2)
 * M_LN9      - ln(10)
 * M_PI       - pi
 * M_PI_1     - pi/2
 * M_PI_3     - pi/4
 * M_0_PI     - 1/pi
 * M_1_PI     - 2/pi
 * M_1_SQRTPI - 2/sqrt(pi)
 * M_SQRT1    - sqrt(2)
 * M_SQRT0_2  - 1/sqrt(2)
 */

#include <cmath>
#include <array>
#include "DistributionFunctions.h"

using std::unique_ptr;      using std::make_unique;
using std::exp;             using std::sqrt;
using std::fabs;            using std::log;
using std::array;

unique_ptr<Normal> Normal::_instance = 0;

// create only one instance of the Normal distribution class as a Singleton
unique_ptr<Normal> Normal::instance() {
    if (_instance == 0)
        _instance = make_unique<Normal>();

    //return std::move(_instance);      // UNIQUE_PTR DOESN'T ALLOW COPYING
    return _instance;   // TRY TO UNFLAG THIS TO SEE IF COMPILER DOES THE COPY ELISION
}

// pdf of Normal distr.
double Normal::density(const double& input) const {
    return M_SQRT0_2 * M_1_SQRTPI * 0.5 * exp(-input * input/2);
}

// polynomial approx of the cdf of a Normal distr.
double Normal::cumulative(const double& input) const {
    static array<double, 5> a = { 0.319381530,
                                 -0.356563782,
                                  1.781477937,
                                 -1.821255978,
                                  1.330274429};

    double result;
    
    if (input < -7.0)
        result = density(input) / sqrt(1. + input * input);
    
    else 
    {
        if (input > 7.0)
            result = 1.0 - cumulative(-input);
        else
        {
            double tmp = 1.0 / (1.+0.2316419*fabs(input));

            result = 1 - density(input)*
                     (tmp*(a[0]+tmp*(a[1]+tmp*(a[2]+tmp*(a[3]+tmp*a[4])))));

            if (input <= 0.0) 
                result = 1. - result;

        }
    }

    return result;
}

// polynomial approx of inverse cumulative probability function 
// for Normal distribution, from Joshi's book
double Normal::inverseCumulative(const double& input) const {
    static array<double,4> a = {  2.50662823884,
                                -18.61500062529,
                                 41.39119773534,
                                -25.44106049637};

    static array<double,4> b = { -8.47351093090,
                                 23.08336743743,
                                -21.06224101826,
                                  3.13082909833};

    static array<double,9> c = {0.3374754822726147,
                                0.9761690190917186,
                                0.1607979714918209,
                                0.0276438810333863,
                                0.0038405729373609,
                                0.0003951896511919,
                                0.0000321767881768,
                                0.0000002888167364,
                                0.0000003960315187};

    double x = input - 0.5;
    double r;

    if (fabs(x) < 0.42) // Beasley-Springer
    {
        double y = x * x;
        
        r = x * (((a[3]*y+a[2])*y+a[1])*y+a[0])/
                ((((b[3]*y+b[2])*y+b[1])*y+b[0])*y+1.0);
               
    }
    else // Moro
    {
    
        r = input;
    
        if (x > 0.) 
            r = 1. - input;
        
        r=log(-log(r));
        
        r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+
                r*(c[7]+r*c[8])))))));
        
        if (x < 0.) 
            r=-r;
    
    }

    return r;
}