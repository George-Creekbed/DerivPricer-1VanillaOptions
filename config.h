#ifndef CONFIG_H
#define CONFIG_H

using Time = double;
using Money = double;

enum class NewtonCotes::Formula;

constexpr NewtonCotes::Formula integration_choice = NewtonCotes::Formula::Simpsons;
constexpr int num_int_intervals = 1000;

#endif // CONFIG_H