#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/riskneutralpricing/heston/heston.h"

namespace py = pybind11;

void PyHeston(py::module &m) {

    py::class_<Heston>(m, "Heston")
        .def(py::init<const double, const double, const double, const double, const double>(),
             py::arg("kappa"), py::arg("theta"), py::arg("eta"), py::arg("rho"), py::arg("v0"))
        .def(py::init<const double, const double, const double, const double, const double, int>(),
             py::arg("kappa"), py::arg("theta"), py::arg("eta"), py::arg("rho"), py::arg("v0"), py::arg("quad_points"))
        .def("get_kappa", &Heston::get_kappa)
        .def("get_theta", &Heston::get_theta)
        .def("get_eta", &Heston::get_eta)
        .def("get_rho", &Heston::get_rho)
        .def("get_v0", &Heston::get_v0)
        .def("get_quadrature", &Heston::get_quadrature)
        .def("get_characteristic_function", &Heston::get_characteristic_function, py::arg("u"), py::arg("T"))
        .def("is_feller_condition_satisfied", &Heston::is_feller_condition_satisfied)
        .def("get_lewis_normalized_price", &Heston::get_lewis_normalized_price, py::arg("x"), py::arg("T"), py::arg("is_call"))
        .def("get_lewis_price", &Heston::get_lewis_price, py::arg("F"), py::arg("K"), py::arg("T"), py::arg("r"), py::arg("is_call"))
        .def("set_kappa", &Heston::set_kappa, py::arg("kappa"))
        .def("set_theta", &Heston::set_theta, py::arg("theta"))
        .def("set_eta", &Heston::set_eta, py::arg("eta"))
        .def("set_rho", &Heston::set_rho, py::arg("rho"))
        .def("set_v0", &Heston::set_v0, py::arg("v0"))
        .def("set_quadrature", &Heston::set_quadrature, py::arg("points"));

    py::enum_<HestonDiscretizationMethod>(m, "HestonDiscretizationMethod")
        .value("EULER", HestonDiscretizationMethod::EULER)
        .value("MILSTEIN", HestonDiscretizationMethod::MILSTEIN)
        .export_values();

    py::class_<HestonSimulation>(m, "HestonSimulation")
        .def(py::init<const double, const double, const double, const double, const double, const double, const double, const double, const int, const int>(),
             py::arg("S"), py::arg("mu"), py::arg("T"), py::arg("kappa"), py::arg("theta"), py::arg("eta"), py::arg("rho"), py::arg("v0"),
             py::arg("number_iterations"), py::arg("sample_size"))
        .def(py::init<const double, const double, const double, const double, const double, const double, const double, const double, HestonDiscretizationMethod, const int, const int>(),
             py::arg("S"), py::arg("mu"), py::arg("T"), py::arg("kappa"), py::arg("theta"), py::arg("eta"), py::arg("rho"), py::arg("v0"),
             py::arg("discretization_method"), py::arg("number_iterations"), py::arg("sample_size"))
        .def("get_S", &HestonSimulation::get_S)
        .def("get_mu", &HestonSimulation::get_mu)
        .def("get_T", &HestonSimulation::get_T)
        .def("get_kappa", &HestonSimulation::get_kappa)
        .def("get_theta", &HestonSimulation::get_theta)
        .def("get_eta", &HestonSimulation::get_eta)
        .def("get_rho", &HestonSimulation::get_rho)
        .def("get_v0", &HestonSimulation::get_v0)
        .def("is_updated", &HestonSimulation::is_updated)
        .def("get_discretization_method", &HestonSimulation::get_discretization_method)
        .def("get_simulated_prices", &HestonSimulation::get_simulated_prices_vector)
        .def("get_simulation_engine", &HestonSimulation::get_simulation_engine)
        .def("get_payoff_engine", &HestonSimulation::get_payoff_engine)
        .def("get_time_taken_for_simulation", &HestonSimulation::get_time_taken_for_simulation)
        .def("get_variance_from_truncation", &HestonSimulation::get_variance_from_truncation, py::arg("v"))
        .def("get_heston", &HestonSimulation::get_heston)
        .def("set_S", &HestonSimulation::set_S, py::arg("S"))
        .def("set_mu", &HestonSimulation::set_mu, py::arg("mu"))
        .def("set_T", &HestonSimulation::set_T, py::arg("T"))
        .def("set_kappa", &HestonSimulation::set_kappa, py::arg("kappa"))
        .def("set_theta", &HestonSimulation::set_theta, py::arg("theta"))
        .def("set_eta", &HestonSimulation::set_eta, py::arg("eta"))
        .def("set_rho", &HestonSimulation::set_rho, py::arg("rho"))
        .def("set_v0", &HestonSimulation::set_v0, py::arg("v0"))
        .def("set_sample_size", &HestonSimulation::set_sample_size, py::arg("sample_size"))
        .def("set_number_iterations", &HestonSimulation::set_number_iterations, py::arg("number_iterations"))
        .def("set_full_truncation", &HestonSimulation::set_full_truncation, py::arg("full_truncation"))
        .def("set_discretization_method", &HestonSimulation::set_discretization_method, py::arg("method"));
}