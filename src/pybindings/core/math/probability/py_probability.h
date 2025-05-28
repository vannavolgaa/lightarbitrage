#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "py_stats.h"
#include "py_distribution.h"
#include "py_simulation.h"

namespace py = pybind11;

void PyProbability(py::module_& m) {

    py::module stats_module = probability_module.def_submodule("stats", "Module for statistical tools");
    py::module simulation_module = probability_module.def_submodule("simulation", "Module for Monte Carlo simulations");
    py::module distribution_module = probability_module.def_submodule("distribution", "Module for probability distributions");

    PyDistribution(distribution_module);
    PySimulation(simulation_module);
    PyStats(stats_module);
}