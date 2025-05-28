#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_errors_frameworks.h"
#include "py_errors_math.h"


namespace py = pybind11;

void PyErrorsCore(py::module_& m)
{
    py::module emath = m.def_submodule("math", "Math errors");
    py::module eframeworks = m.def_submodule("frameworks", "Frameworks errors");
    PyErrorsFrameworks(eframeworks);
    PyErrorsMath(emath);
    
}