#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/parametrization/svi/svi.h"

namespace py = pybind11;

void PySVI(py::module_& m)
{
    py::class_<SSVI, std::shared_ptr<SSVI>>(m, "SSVI")
        .def(py::init<double, double, double>(), py::arg("rho"), py::arg("nu"), py::arg("gamma"))
        .def("set_rho", &SSVI::set_rho, py::arg("rho"))
        .def("set_nu", &SSVI::set_nu, py::arg("nu"))
        .def("set_gamma", &SSVI::set_gamma, py::arg("gamma"))
        .def("get_rho", &SSVI::get_rho)
        .def("get_nu", &SSVI::get_nu)
        .def("get_gamma", &SSVI::get_gamma)
        .def("get_svi", &SSVI::get_svi, py::arg("atm_total_variance"), py::arg("t"))
        .def("check_butterfly_arbitrage", &SSVI::check_butterfly_arbitrage, py::arg("atm_total_variance"))
        .def("check_calendar_spread_arbitrage", &SSVI::check_calendar_spread_arbitrage, py::arg("atm_total_variance"))
        .def("get_total_variance", &SSVI::get_total_variance, py::arg("k"), py::arg("atm_total_variance"))
        .def("get_implied_variance", &SSVI::get_implied_variance, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"))
        .def("get_implied_volatility", &SSVI::get_implied_volatility, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"))
        .def("get_risk_neutral_density", &SSVI::get_risk_neutral_density, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"))
        .def("get_local_volatility", &SSVI::get_local_volatility, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"))
        .def("get_atm_volatility_skew", &SSVI::get_atm_volatility_skew, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"))
        .def("get_undiscounted_normalized_black_price", &SSVI::get_undiscounted_normalized_black_price, py::arg("k"), py::arg("atm_total_variance"), py::arg("t"), py::arg("is_call"), "Compute the undiscounted normalized Black price");

    py::class_<SVI, std::shared_ptr<SVI>>(m, "SVI")
        .def(py::init<double, double, double, double, double, double>(), 
             py::arg("vt"), py::arg("ut"), py::arg("ct"), py::arg("pt"), py::arg("vmt"), py::arg("t"))
        .def("get_vt", &SVI::get_vt)
        .def("get_ut", &SVI::get_ut)
        .def("get_ct", &SVI::get_ct)
        .def("get_pt", &SVI::get_pt)
        .def("get_vmt", &SVI::get_vmt)
        .def("get_T", &SVI::get_T)
        .def("set_vt", &SVI::set_vt, py::arg("vt"))
        .def("set_ut", &SVI::set_ut, py::arg("ut"))
        .def("set_ct", &SVI::set_ct, py::arg("ct"))
        .def("set_pt", &SVI::set_pt, py::arg("pt"))
        .def("set_vmt", &SVI::set_vmt, py::arg("vmt"))
        .def("set_T", &SVI::set_T, py::arg("T"))
        .def("get_ssvi", &SVI::get_ssvi)
        .def("check_butterfly_arbitrage", &SVI::check_butterfly_arbitrage)
        .def("check_calendar_spread_arbitrage", &SVI::check_calendar_spread_arbitrage, py::arg("slice"), py::arg("grid_size"), py::arg("epsilon"))
        .def("get_total_variance", &SVI::get_total_variance, py::arg("k"))
        .def("get_atm_total_variance", &SVI::get_atm_total_variance)
        .def("get_minimum_total_variance", &SVI::get_minimum_total_variance)
        .def("get_implied_variance", &SVI::get_implied_variance, py::arg("k"))
        .def("get_implied_volatility", &SVI::get_implied_volatility, py::arg("k"))
        .def("get_risk_neutral_density", &SVI::get_risk_neutral_density, py::arg("k"))
        .def("get_local_variance", &SVI::get_local_variance, py::arg("k"))
        .def("get_local_volatility", &SVI::get_local_volatility, py::arg("k"))
        .def("get_asymptotic_value", &SVI::get_asymptotic_value, py::arg("left"), py::arg("epsilon"))
        .def("get_undiscounted_normalized_black_price", &SVI::get_undiscounted_normalized_black_price, py::arg("k"), py::arg("t"), py::arg("is_call"), "Compute the undiscounted normalized Black price");

    py::class_<ReducedSVI, SVI, std::shared_ptr<ReducedSVI>>(m, "ReducedSVI")
        .def(py::init<double, double, double, double>(), py::arg("vt"), py::arg("eta"), py::arg("rho"), py::arg("t"))
        .def("get_eta", &ReducedSVI::get_eta)
        .def("get_rho", &ReducedSVI::get_rho)
        .def("set_eta", &ReducedSVI::set_eta, py::arg("eta"))
        .def("set_rho", &ReducedSVI::set_rho, py::arg("rho"));
}
