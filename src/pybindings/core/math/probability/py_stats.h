#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "../../../../../src/core/math/probability/stats.h"

namespace py = pybind11;

void PyStats(py::module_& m) 
{
    py::class_<CorrelationMatrix, std::shared_ptr<CorrelationMatrix>>(m, "CorrelationMatrix")
        .def(py::init<>())
        .def(py::init<double>())
        .def(py::init<const std::vector<double>&>())
        .def(py::init<const Eigen::MatrixXd&>())
        .def("get_matrix", &CorrelationMatrix::get_matrix)
        .def("get_cholesky_decomposition", &CorrelationMatrix::get_cholesky_decomposition)
        .def("get_dimension", &CorrelationMatrix::get_dimension);

    // Add more bindings for other classes as needed
}