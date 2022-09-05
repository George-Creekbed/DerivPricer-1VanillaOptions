#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <utility>      // std::move
#include <memory>       // std::unique_ptr
#include <functional>   // std::function
#include <algorithm>    // std::copy
#include <numeric>      // std::accumulate
#include "config.h"
#include "NewtonCotesFormulas.h"

template<class T>
class IntegrationStrategy {
public:
    virtual T integrate(const Time, const Time) const = 0;
    virtual ~IntegrationStrategy(){}
};

template<class T>
class AnalyticIntegration: public IntegrationStrategy {
public:
    AnalyticIntegration(const std::function<T (const Time&)>&);
    T integrate(const Time, const Time) const;
private:
    std::function<T (const Time&)> already_integrated_function; // lambda of integral(\sigma(S(t))), for instance
};

template<class T>
class NumericIntegration: public IntegrationStrategy {
public:
    NumericIntegration(const std::function<T (const Time&)>&);
    T integrate(const Time, const Time) const;
private:
    std::function<T (const Time&)> not_yet_integrated_function; // lambda of integral(\sigma(S(t))), for instance
};

template<class T, template<class T> class container>
class DiscreteDataIntegration: public IntegrationStrategy {
public:
    DiscreteDataIntegration(typename container<T>::const_iterator&, typename container<T>::const_iterator&);
    T integrate(const Time, const Time) const;
private:
    container<T> data_points;
};

// move-only class. Copy by creating new pointer with clone()
template<class T> 
class Parameter {
public:
    Parameter(){}
    Parameter(const T, const IntegrationStrategy<T>&); //, const IntegrationType&);
    virtual ~Parameter(){}

    Parameter<T>(const Parameter<T>&) = delete;
    Parameter<T>& operator=(const Parameter<T>&) = delete;
    Parameter<T>(Parameter<T>&&) noexcept; // = default;
    Parameter<T>& operator=(Parameter<T>&&) noexcept; // = default;

    virtual std::unique_ptr<T> clone() const;
    virtual T integrate(const Time, const Time) const;
    T mean(const Time, const Time) const;
    
protected:
    T object;
    IntegrationStrategy<T> strategy;
};

template<class T> AnalyticIntegration<T>::AnalyticIntegration(const std::function<T (const Time&)>& function_):
    already_integrated_function(function_) {}

template<class T> T AnalyticIntegration<T>::integrate(const Time t1, const Time t2) const {
    return already_integrated_function(t2) - already_integrated_function(t1);
}

template<class T> NumericIntegration<T>::NumericIntegration(const std::function<T (const Time&)>& function_):
    not_yet_integrated_function(function_) {}

template<class T> T NumericIntegration<T>::integrate(const Time t1, const Time t2) const {
    return NewtonCotes::createQuadrature(integration_choice)(t1, t2, num_int_intervals, not_yet_integrated_function);
}

template<class T, template<class T> class container> 
DiscreteDataIntegration<T, container>::DiscreteDataIntegration(typename container<T>::const_iterator& begin, typename container<T>::const_iterator& end) {
    std::copy(begin, end, data_points.begin());
}
    
template<class T, template<class T> class container> 
T DiscreteDataIntegration<T, container>::integrate(const Time t1, const Time t2) const {
    double weight = (t2 - t1) / data_points.size();
    return weight * std::accumulate(data_points.begin(), data_points.end(), 0.0);
}


template<class T> Parameter<T>::Parameter(const T original, const IntegrationStrategy<T>& strategy_): 
    object(original), strategy(strategy_) {}

/* Copy disallowed
template<class T> Parameter<T>::Parameter(const Parameter<T>& original) {
    object = original.object;
}

template<class T> Parameter<T>& Parameter<T>::operator=(const Parameter<T>& original) {
    if (&original != this)
        object = original.object;
    
    return *this;
} */

template<class T> Parameter<T>::Parameter(Parameter<T>&& original) noexcept {
        object = std::move(original.object); // std::move is noexcept guaranteed
        strategy = std::move(original.strategy);
}

template<class T> Parameter<T>& Parameter<T>::operator=(Parameter<T>&& original) noexcept {
    if (&original != this) {
        object = std::move(original.object);
        strategy = std::move(original.strategy);
    }
    return *this;
}

template<class T> std::unique_ptr<T> Parameter<T>::clone() const {
    // deep copy
    // just an alias for make_unique at this point, ie just use make_unique directly
    return std::make_unique<T>(*this);
}

template<class T> T Parameter<T>::integrate(const Time t1, const Time t2) const {
    return strategy.integrate(t1, t2);
}

template<class T> T Parameter<T>::mean(const Time t1, const Time t2) const {
    return integrate(t1, t2) / (t2 - t1);
}

//     friend template<class T> 
//     T integral(const Time, const Time,
//                     std::function<T (const Time, const Time, const T&)>);

//     T mean(const Time, const Time) const;
//     T rootMeanSquare(const Time t_in, const Time t_fin) const;

// private:
//     T object;
//     IntegrationType int_choice = IntegrationType::analytic; 
// };

// template<class T>
// T integral(const Time t_in, const Time t_fin,
//                 std::function<T (const Time, const Time, const T&)> int_fun) {
//     // variables integration_choice and num_int_intervals are in "config.h"
//     if (int_choice == IntegrationType::analytic)
//         return int_fun(t_in, t_fin, object);
//     else if (int_choice == IntegrationType::numeric) {
//         auto num_int = NewtonCotes::createQuadrature(integration_choice);
//         return num_int(t_in, t_fin, num_int_intervals, int_fun);
//     } else if (int_choice == IntegrationType::datapoints)
//         return integral(data.begin(), data.end());
//     else
//         throw std::logic_error("IntegrationType should be either analytic or numeric.\n");

//     return 0;
// }

// template<class T>
// T integral(const Time t_in, const Time t_fin, 

// template<class T>
// T Parameter<T>::mean(const Time t_in, const Time t_fin) const {
//     return integral(t_in, t_fin )
//}

#endif // PARAMETERS_H