#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_svi.h"
#include "py_nss.h"

namespace py = pybind11;

void PyParametrization(py::module_& m)
{
    py::module svi = m.def_submodule("svi", "Definition of stochastic voatility inspired model from Gatheral.");
    py::module nss = m.def_submodule("nelsonsiegelsvensson", "Definition of Nelson Siegel and Nelson Siegel Svensson models");
    PyNSS(nss);
    PySVI(svi);
}

