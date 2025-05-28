#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "../../../../../src/core/math/probability/simulation.h"

namespace py = pybind11;

void PySimulation(py::module_& m) 
{

    py::class_<MonteCarloEngine, std::shared_ptr<MonteCarloEngine>>(m, "MonteCarloEngine")
        .def(py::init<int, int, std::shared_ptr<ProbabilityDistribution>, bool>(), 
             py::arg("number_iterations"), py::arg("sample_size"), py::arg("distribution"), py::arg("antithetic"))
        .def(py::init<int, std::shared_ptr<ProbabilityDistribution>, bool>(), 
             py::arg("number_iterations"), py::arg("distribution"), py::arg("antithetic"))
        .def("get_samples", &MonteCarloEngine::get_samples)
        .def("get_sample_means", &MonteCarloEngine::get_sample_means)
        .def("get_sample_variances", &MonteCarloEngine::get_sample_variances)
        .def("get_sample_size", &MonteCarloEngine::get_sample_size)
        .def("get_number_iterations", &MonteCarloEngine::get_number_iterations)
        .def("get_distribution", &MonteCarloEngine::get_distribution)
        .def("is_updated", &MonteCarloEngine::is_updated)
        .def("is_antithetic", &MonteCarloEngine::is_antithetic)
        .def("get_time_taken_for_samples_update", &MonteCarloEngine::get_time_taken_for_samples_update)
        .def("set_sample_dimension", &MonteCarloEngine::set_sample_dimension)
        .def("set_distribution", &MonteCarloEngine::set_distribution)
        .def("set_antithetic", &MonteCarloEngine::set_antithetic)
        .def("update_samples", &MonteCarloEngine::update_samples);

    py::class_<MultivariateMonteCarloEngine, std::shared_ptr<MultivariateMonteCarloEngine>>(m, "MultivariateMonteCarloEngine")
        .def(py::init<const CorrelationMatrix&, int, int, std::shared_ptr<ProbabilityDistribution>>(),
             py::arg("correlation_matrix"), py::arg("number_iterations"), py::arg("sample_size"), py::arg("distribution"))
        .def("get_samples", &MultivariateMonteCarloEngine::get_samples)
        .def("get_sample_means", &MultivariateMonteCarloEngine::get_sample_means)
        .def("get_sample_variances", &MultivariateMonteCarloEngine::get_sample_variances)
        .def("get_sample_size", &MultivariateMonteCarloEngine::get_sample_size)
        .def("get_number_iterations", &MultivariateMonteCarloEngine::get_number_iterations)
        .def("get_correlation_matrix", &MultivariateMonteCarloEngine::get_correlation_matrix)
        .def("get_distribution", &MultivariateMonteCarloEngine::get_distribution)
        .def("is_updated", &MultivariateMonteCarloEngine::is_updated)
        .def("get_time_taken_for_samples_update", &MultivariateMonteCarloEngine::get_time_taken_for_samples_update)
        .def("set_sample_dimension", &MultivariateMonteCarloEngine::set_sample_dimension)
        .def("set_correlation_matrix", &MultivariateMonteCarloEngine::set_correlation_matrix)
        .def("set_distribution", &MultivariateMonteCarloEngine::set_distribution)
        .def("update_samples", &MultivariateMonteCarloEngine::update_samples);
}