#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/riskneutralpricing/blackscholes/bs.h"

namespace py = pybind11;

void PyBlackScholes(py::module &m) {

    py::class_<BlackScholes>(m, "BlackScholes")
        .def(py::init<const double, const double, const double, const double, const double, const double, const bool, const bool>(),
             py::arg("S"), py::arg("K"), py::arg("T"), py::arg("r"), py::arg("q"), py::arg("sigma"), py::arg("is_call"), py::arg("is_future"))
        .def("get_S", &BlackScholes::get_S)
        .def("get_K", &BlackScholes::get_K)
        .def("get_T", &BlackScholes::get_T)
        .def("get_r", &BlackScholes::get_r)
        .def("get_q", &BlackScholes::get_q)
        .def("get_sigma", &BlackScholes::get_sigma)
        .def("get_is_call", &BlackScholes::get_is_call)
        .def("get_is_future", &BlackScholes::get_is_future)
        .def("set_S", &BlackScholes::set_S, py::arg("S"))
        .def("set_K", &BlackScholes::set_K, py::arg("K"))
        .def("set_T", &BlackScholes::set_T, py::arg("T"))
        .def("set_r", &BlackScholes::set_r, py::arg("r"))
        .def("set_q", &BlackScholes::set_q, py::arg("q"))
        .def("set_sigma", &BlackScholes::set_sigma, py::arg("sigma"))
        .def("set_is_call", &BlackScholes::set_is_call, py::arg("is_call"))
        .def("set_is_future", &BlackScholes::set_is_future, py::arg("is_future"))
        .def("get_price", &BlackScholes::get_price)
        .def("get_delta", &BlackScholes::get_delta)
        .def("get_gamma", &BlackScholes::get_gamma)
        .def("get_vega", &BlackScholes::get_vega)
        .def("get_theta", &BlackScholes::get_theta)
        .def("get_dual_delta", &BlackScholes::get_dual_delta)
        .def("get_dual_gamma", &BlackScholes::get_dual_gamma);

    py::class_<BaroneAdesiWhaley>(m, "BaroneAdesiWhaley")
        .def(py::init<const double, const double, const double, const double, const double, const double, const bool, const bool>(),
             py::arg("S"), py::arg("K"), py::arg("T"), py::arg("r"), py::arg("q"), py::arg("sigma"), py::arg("is_call"), py::arg("is_future"))
        .def("get_S", &BaroneAdesiWhaley::get_S)
        .def("get_K", &BaroneAdesiWhaley::get_K)
        .def("get_T", &BaroneAdesiWhaley::get_T)
        .def("get_r", &BaroneAdesiWhaley::get_r)
        .def("get_q", &BaroneAdesiWhaley::get_q)
        .def("get_sigma", &BaroneAdesiWhaley::get_sigma)
        .def("get_is_call", &BaroneAdesiWhaley::get_is_call)
        .def("get_is_future", &BaroneAdesiWhaley::get_is_future)
        .def("set_S", &BaroneAdesiWhaley::set_S, py::arg("S"))
        .def("set_K", &BaroneAdesiWhaley::set_K, py::arg("K"))
        .def("set_T", &BaroneAdesiWhaley::set_T, py::arg("T"))
        .def("set_r", &BaroneAdesiWhaley::set_r, py::arg("r"))
        .def("set_q", &BaroneAdesiWhaley::set_q, py::arg("q"))
        .def("set_sigma", &BaroneAdesiWhaley::set_sigma, py::arg("sigma"))
        .def("set_is_call", &BaroneAdesiWhaley::set_is_call, py::arg("is_call"))
        .def("set_is_future", &BaroneAdesiWhaley::set_is_future, py::arg("is_future"))
        .def("get_exercise_premium", &BaroneAdesiWhaley::get_exercise_premium)
        .def("get_price", &BaroneAdesiWhaley::get_price);

    py::class_<BlackScholesSimulation>(m, "BlackScholesSimulation")
        .def(py::init<const double, const double, const double, const double, const int, const int, bool>(),
             py::arg("S"), py::arg("mu"), py::arg("sigma"), py::arg("T"), py::arg("number_iterations"), py::arg("sample_size"), py::arg("antithetic"))
        .def("get_S", &BlackScholesSimulation::get_S)
        .def("get_mu", &BlackScholesSimulation::get_mu)
        .def("get_sigma", &BlackScholesSimulation::get_sigma)
        .def("get_T", &BlackScholesSimulation::get_T)
        .def("get_simulated_prices", &BlackScholesSimulation::get_simulated_prices_vector)
        .def("get_time_taken_for_simulation", &BlackScholesSimulation::get_time_taken_for_simulation)
        .def("set_S", &BlackScholesSimulation::set_S, py::arg("S"))
        .def("set_mu", &BlackScholesSimulation::set_mu, py::arg("mu"))
        .def("set_sigma", &BlackScholesSimulation::set_sigma, py::arg("sigma"))
        .def("set_T", &BlackScholesSimulation::set_T, py::arg("T"))
        .def("set_sample_size", &BlackScholesSimulation::set_sample_size, py::arg("sample_size"))
        .def("set_number_iterations", &BlackScholesSimulation::set_number_iterations, py::arg("number_iterations"))
        .def("get_simulation_engine", &BlackScholesSimulation::get_simulation_engine)
        .def("get_payoff_engine", &BlackScholesSimulation::get_payoff_engine)
        .def("is_updated", &BlackScholesSimulation::is_updated);

}