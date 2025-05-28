#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/core/math/erf_cody.h"

namespace py = pybind11;

void PyErfCody(py::module_& m)
{
    m.def("erf_cody", &erf_cody, py::arg("x"), "Accurate implementation of erf(x)");
    m.def("erfc_cody", &erfc_cody, py::arg("x"), "Accurate implementation of erfc(x)");
    m.def("erfcx_cody", &erfcx_cody, py::arg("x"), "Accurate implementation of scaled erfcx(x) = exp(x^2) * erfc(x)");
}

