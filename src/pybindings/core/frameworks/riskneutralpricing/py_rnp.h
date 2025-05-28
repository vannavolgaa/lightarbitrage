#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_blackscholes.h"
#include "py_heston.h"
#include "py_lbr.h"
#include "py_payoff.h"

namespace py = pybind11;

void PyRNP(py::module &m) {

    py::module bs = m.def_submodule("blackscholes", "Black-Scholes framework for risk-neutral pricing.");
    py::module heston = m.def_submodule("heston", "Heston model for risk-neutral pricing.");
    py::module lbr = m.def_submodule("letsberational", "LetsBeRational framework for implied volatility calculation.");
    py::module payoff = m.def_submodule("payoff", "Payoff engine for risk-neutral pricing.");

    PyBlackScholes(bs);
    PyHeston(heston);
    PyLBR(lbr);
    PyPayoff(payoff);
    
}