#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "../../../../src/core/math/optimization/optimization.h"

namespace py = pybind11;

class PyTrampolineOptimizer : public Optimizer {
public:
    using Optimizer::Optimizer;

    std::shared_ptr<OptimizerResult> optimize() const override {
        PYBIND11_OVERRIDE_PURE(std::shared_ptr<OptimizerResult>, Optimizer, optimize);
    }
};

void PyOptimization(py::module_& m) {

    py::class_<OptimizerResult, std::shared_ptr<OptimizerResult>>(m, "OptimizerResult")
        .def_readonly("x", &OptimizerResult::x_)
        .def_readonly("function_value", &OptimizerResult::function_value_)
        .def_readonly("number_iterations", &OptimizerResult::number_iterations_)
        .def_readonly("tolerance_threshold", &OptimizerResult::tolerance_threshold_)
        .def_readonly("maximum_iterations", &OptimizerResult::maximum_iterations_)
        .def_readonly("time_taken", &OptimizerResult::time_taken_);

    py::class_<Optimizer, PyTrampolineOptimizer, std::shared_ptr<Optimizer>>(m, "Optimizer")
        .def("optimize", &Optimizer::optimize);

    py::class_<NewtonRaphson, Optimizer, std::shared_ptr<NewtonRaphson>>(m, "NewtonRaphson")
        .def(py::init<double, std::function<double(double)>, std::function<double(double)>>(),
                py::arg("x0"), py::arg("f"), py::arg("df"))
        .def("set_tolerance_rate", &NewtonRaphson::set_tolerance_rate)
        .def("set_max_iterations", &NewtonRaphson::set_max_iterations)
        .def("optimize", &NewtonRaphson::optimize);

    // Nelder-Mead
    py::class_<NelderMead, Optimizer, std::shared_ptr<NelderMead>>(m, "NelderMead")
        .def(py::init<std::vector<double>, std::function<double(std::vector<double>)>>(),
                py::arg("x0"), py::arg("f"))
        .def("set_initial_simplex_method", &NelderMead::set_initial_simplex_method)
        .def("set_max_iterations", &NelderMead::set_max_iterations)
        .def("set_perturbation_parameter", &NelderMead::set_perturbation_parameter)
        .def("set_reflection_parameter", &NelderMead::set_reflection_parameter)
        .def("set_expansion_parameter", &NelderMead::set_expansion_parameter)
        .def("set_contraction_parameter", &NelderMead::set_contraction_parameter)
        .def("set_shrink_parameter", &NelderMead::set_shrink_parameter)
        .def("set_tolerance_rate", &NelderMead::set_tolerance_rate)
        .def("optimize", &NelderMead::optimize);
}
    