#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/core/math/quadratures.h"

namespace py = pybind11;

void PyQuadratures(py::module_& m)
{
    py::class_<GaussianQuadrature, std::shared_ptr<GaussianQuadrature>>(m, "GaussianQuadrature")
        .def(py::init<int, const std::vector<double>&, const std::vector<double>&>(),
             py::arg("points"), py::arg("roots"), py::arg("weights"))
        .def_readonly("points", &GaussianQuadrature::points)
        .def_readonly("roots", &GaussianQuadrature::roots)
        .def_readonly("weights", &GaussianQuadrature::weights)
        .def("integrate", &GaussianQuadrature::integrate, py::arg("f"),
             "Compute the Gaussian quadrature integral of a function.");

    m.def("get_gauss_laguerre_quadrature", &get_gauss_laguerre_quadrature, py::arg("points"),
          "Return a Gaussian quadrature object for Gauss-Laguerre integration.");
    
}

