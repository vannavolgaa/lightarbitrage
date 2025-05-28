#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../src/pybindings/core/frameworks/py_frameworks.h"
#include "../../../src/pybindings/core/math/py_math.h"

namespace py = pybind11;

void PyCore(py::module_& m)
{
    py::module frameworks = m.def_submodule("frameworks", "Definition of frameworks/model.");
    py::module math = m.def_submodule("math", "Definition of math tools.");

    PyFrameworks(frameworks);
    PyMath(math);

}