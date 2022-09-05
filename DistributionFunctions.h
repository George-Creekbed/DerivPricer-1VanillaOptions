#ifndef DISTRIBUTIONFUNCTIONS_H
#define DISTRIBUTIONFUNCTIONS_H

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