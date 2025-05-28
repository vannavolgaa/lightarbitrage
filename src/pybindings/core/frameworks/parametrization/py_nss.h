#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../../src/core/frameworks/parametrization/nelsonsiegelsvensson/nss.h"

namespace py = pybind11;

void PyNSS(py::module_& m)
{
    py::class_<NelsonSiegel, std::shared_ptr<NelsonSiegel>>(m, "NelsonSiegel")
        .def(py::init<double, double, double, double>(), 
             py::arg("beta0"), py::arg("beta1"), py::arg("beta2"), py::arg("tau"))
        .def("get_beta0", &NelsonSiegel::get_beta0)
        .def("get_beta1", &NelsonSiegel::get_beta1)
        .def("get_beta2", &NelsonSiegel::get_beta2)
        .def("get_tau", &NelsonSiegel::get_tau)
        .def("set_beta0", &NelsonSiegel::set_beta0, py::arg("beta0"))
        .def("set_beta1", &NelsonSiegel::set_beta1, py::arg("beta1"))
        .def("set_beta2", &NelsonSiegel::set_beta2, py::arg("beta2"))
        .def("set_tau", &NelsonSiegel::set_tau, py::arg("tau"))
        .def("get_rate", &NelsonSiegel::get_rate, py::arg("t"))
        .def("get_forward_rate", &NelsonSiegel::get_forward_rate, py::arg("t"));

    py::class_<NelsonSiegelSvensson, NelsonSiegel, std::shared_ptr<NelsonSiegelSvensson>>(m, "NelsonSiegelSvensson")
        .def(py::init<double, double, double, double, double, double>(), 
             py::arg("beta0"), py::arg("beta1"), py::arg("beta2"), py::arg("beta3"), py::arg("tau1"), py::arg("tau2"))
        .def("get_beta3", &NelsonSiegelSvensson::get_beta3)
        .def("get_tau2", &NelsonSiegelSvensson::get_tau2)
        .def("set_beta3", &NelsonSiegelSvensson::set_beta3, py::arg("beta3"))
        .def("set_tau2", &NelsonSiegelSvensson::set_tau2, py::arg("tau2"))
        .def("get_rate", &NelsonSiegelSvensson::get_rate, py::arg("t"))
        .def("get_forward_rate", &NelsonSiegelSvensson::get_forward_rate, py::arg("t"));
    
}
