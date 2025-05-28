#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/riskneutralpricing/payoff/montecarlo.h"

namespace py = pybind11;

void PyPayoff(py::module &m) 
{
    py::class_<MonteCarloPayoffEngine>(m, "MonteCarloPayoffEngine")
        .def(py::init<std::vector<std::vector<double>>, double>(),
             py::arg("prices"), py::arg("T"),
             "Initialize with a vector of vectors representing prices and time to maturity.")
        .def("european_vanilla_price", &MonteCarloPayoffEngine::european_vanilla_price,
             py::arg("r"), py::arg("K"), py::arg("is_call"),
             "Calculate the European vanilla option price.");
}
    
