#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/core/math/interpolation2D.h"

namespace py = pybind11;

class TrampolineInterpolation2D : public Interpolation2D
{
    public: 
        using Interpolation2D::Interpolation2D; 
        double implemented_evaluate(double x_value) const override { PYBIND11_OVERRIDE_PURE(double, Interpolation2D, implemented_evaluate, x_value);}
};

void PyInterpolation2D(py::module_& m)
{
    py::class_<Interpolation2D, TrampolineInterpolation2D, std::shared_ptr<Interpolation2D>>(m, "Interpolation2D")
        .def("set_values", &Interpolation2D::set_values, py::arg("mapped_x_y_"))
        .def("evaluate", &Interpolation2D::evaluate, py::arg("x_value"))
        .def("get_x_min", &Interpolation2D::get_x_min)
        .def("get_x_max", &Interpolation2D::get_x_max)
        .def("get_x_vector", &Interpolation2D::get_x_vector)
        .def("get_y_vector", &Interpolation2D::get_y_vector);

    py::class_<LinearInterpolation2D, Interpolation2D, std::shared_ptr<LinearInterpolation2D>>(m, "LinearInterpolation2D")
        .def(py::init<std::map<double, double>>(), py::arg("mapped_x_y_"));

    py::class_<CubicSpline2D, Interpolation2D, std::shared_ptr<CubicSpline2D>>(m, "CubicSpline2D")
        .def(py::init<std::map<double, double>>(), py::arg("mapped_x_y_"));
}
