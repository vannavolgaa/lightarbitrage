#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "../../../../../src/core/math/probability/probability.h"

namespace py = pybind11;

// Trampoline for ProbabilityDistribution
class PyProbabilityDistribution : public ProbabilityDistribution {
    public:
        using ProbabilityDistribution::ProbabilityDistribution;
        double inv_cdf(double p) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, inv_cdf, p); }
        double cdf(double x) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, cdf, x); }
        double mgf(double t) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, mgf, t); }
        std::complex<double> cf(double t) const override { PYBIND11_OVERRIDE_PURE(std::complex<double>, ProbabilityDistribution, cf, t); }
        double r() const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, r); }
    };
    
    // Trampoline for ContinuousProbabilityDistribution
class PyContinuousProbabilityDistribution : public ContinuousProbabilityDistribution {
    public:
        using ContinuousProbabilityDistribution::ContinuousProbabilityDistribution;
        double inv_cdf(double p) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, inv_cdf, p); }
        double cdf(double x) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, cdf, x); }
        double mgf(double t) const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, mgf, t); }
        std::complex<double> cf(double t) const override { PYBIND11_OVERRIDE_PURE(std::complex<double>, ProbabilityDistribution, cf, t); }
        double r() const override { PYBIND11_OVERRIDE_PURE(double, ProbabilityDistribution, r); }
        double pdf(double x) const override { PYBIND11_OVERRIDE_PURE(double, ContinuousProbabilityDistribution, pdf, x); }
    };

void PyDistribution(py::module_& m) {

    py::class_<ProbabilityDistribution, PyProbabilityDistribution, std::shared_ptr<ProbabilityDistribution>>(m, "ProbabilityDistribution")
        .def("inv_cdf", &ProbabilityDistribution::inv_cdf)
        .def("cdf", &ProbabilityDistribution::cdf)
        .def("mgf", &ProbabilityDistribution::mgf)
        .def("cf", &ProbabilityDistribution::cf)
        .def("r", &ProbabilityDistribution::r);

    py::class_<ContinuousProbabilityDistribution, ProbabilityDistribution, PyContinuousProbabilityDistribution, std::shared_ptr<ContinuousProbabilityDistribution>>(m, "ContinuousProbabilityDistribution")
        .def("pdf", &ContinuousProbabilityDistribution::pdf);

    // Bind NormalDistribution as a concrete class
    py::class_<NormalDistribution, ContinuousProbabilityDistribution, std::shared_ptr<NormalDistribution>>(m, "NormalDistribution")
        .def(py::init<>()) // Default constructor
        .def(py::init<double, double>(), py::arg("mu"), py::arg("sigma")) // Parameterized constructor
        .def("inv_cdf", &NormalDistribution::inv_cdf)
        .def("cdf", &NormalDistribution::cdf)
        .def("mgf", &NormalDistribution::mgf)
        .def("cf", &NormalDistribution::cf)
        .def("pdf", &NormalDistribution::pdf)
        .def("r", &NormalDistribution::r)
        .def("get_mu", &NormalDistribution::get_mu)
        .def("set_mu", &NormalDistribution::set_mu)
        .def("get_sigma", &NormalDistribution::get_sigma)
        .def("set_sigma", &NormalDistribution::set_sigma);
    
    py::class_<UniformDistribution, ContinuousProbabilityDistribution, std::shared_ptr<UniformDistribution>>(m, "UniformDistribution")
        .def(py::init<double, double>(), py::arg("left_bound"), py::arg("right_bound"))
        .def(py::init<>()) // Default constructor
        .def("inv_cdf", &UniformDistribution::inv_cdf)
        .def("cdf", &UniformDistribution::cdf)
        .def("mgf", &UniformDistribution::mgf)
        .def("cf", &UniformDistribution::cf)
        .def("pdf", &UniformDistribution::pdf)
        .def("r", &UniformDistribution::r)
        .def("get_left_bound", &UniformDistribution::get_left_bound)
        .def("get_right_bound", &UniformDistribution::get_right_bound);
        .def("set_bounds", &UniformDistribution::set_bounds);

}
