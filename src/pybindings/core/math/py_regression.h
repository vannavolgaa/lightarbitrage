#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/core/math/regression.h"

namespace py = pybind11;

void PyRegression(py::module_& m)
{
    py::class_<RegressionResult>(m, "RegressionResult")
        .def_readonly("intercept", &RegressionResult::intercept)
        .def_readonly("coefficients", &RegressionResult::coefficients)
        .def_readonly("residuals", &RegressionResult::residuals)
        .def_readonly("r_squared", &RegressionResult::r_squared)
        .def_readonly("rmse", &RegressionResult::rmse)
        .def_readonly("mae", &RegressionResult::mae)
        .def("get_predicted_value", &RegressionResult::get_predicted_value, py::arg("x"),
             "Get the predicted value for a given input vector.");

    // Bind the compute_linear_least_square function
    m.def("compute_linear_least_square", &compute_linear_least_square, py::arg("Y"), py::arg("X"), py::arg("use_intercept"),"Compute linear least squares regression.");
}