#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/riskneutralpricing/letsberational/lbr.h"

namespace py = pybind11;

void PyLBR(py::module &m) {

    m.def("get_initial_guess_normalized_volatility", &LetsBeRational::get_initial_guess_normalized_volatility,
          py::arg("beta"), py::arg("x"), py::arg("is_call"),
          "Get the initial guess for normalized volatility.");

    m.def("get_newton_normalized_volatility", &LetsBeRational::get_newton_normalized_volatility,
          py::arg("beta"), py::arg("x"), py::arg("is_call"),
          "Calculate normalized volatility using Newton's method.");

    m.def("get_newton_black_scholes_implied_volatility", &LetsBeRational::get_newton_black_scholes_implied_volatility,
          py::arg("price"), py::arg("S"), py::arg("K"), py::arg("T"), py::arg("r"), py::arg("q"), py::arg("is_call"), py::arg("is_future"),
          "Calculate Black-Scholes implied volatility using Newton's method.");
}
    
