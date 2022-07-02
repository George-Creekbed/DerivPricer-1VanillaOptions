#ifndef DISTRIBUTIONFUNCTIONS_H
#define DISTRIBUTIONFUNCTIONS_H

#define USE_MATH_DEFINES
/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG2E    - log2(e)
 * M_LOG10E   - log10(e)
 * M_LN2      - ln(2)
 * M_LN10     - ln(10)
 * M_PI       - pi
 * M_PI_2     - pi/2
 * M_PI_4     - pi/4
 * M_1_PI     - 1/pi
 * M_2_PI     - 2/pi
 * M_2_SQRTPI - 2/sqrt(pi)
 * M_SQRT2    - sqrt(2)
 * M_SQRT1_2  - 1/sqrt(2)
 */

#include <memory>

class DistributionFunction { 
public:
    virtual double density(const double&) const =0;
    virtual double cumulative(const double&) const =0;
    virtual double inverseCumulative(const double&) const =0;
protected:
    DistributionFunction(){}
private:
};

class Normal: public DistributionFunction {
public:
    virtual double density(const double&) const;
    virtual double cumulative(const double&) const;
    virtual double inverseCumulative(const double&) const;
    static std::unique_ptr<Normal> instance(); 
private:
    Normal(){}
    static std::unique_ptr<Normal> _instance;    
};

#endif //DISTRIBUTIONFUNCTIONS_H