#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/pybindings/core/frameworks/parametrization/py_parametrization.h"
#include "../../../../src/pybindings/core/frameworks/riskneutralpricing/py_rnp.h"

namespace py = pybind11;

void PyFrameworks(py::module_& m)
{
    py::module parametrization = m.def_submodule("parametrization", "Definition of parametrization frameworks");
    py::module rnp = m.def_submodule("riskneutralpricing", "Definition of models for risk neutral pricing.");

    PyParametrization(parametrization);
    PyRNP(rnp);

}