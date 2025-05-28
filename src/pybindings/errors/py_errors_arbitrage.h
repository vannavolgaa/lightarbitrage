#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../src/pybindings/errors/core/py_errors_core.h"
#include "../../../src/errors/base.h"

namespace py = pybind11;

void PyErrorsArbitrage(py::module_& m)
{
    py::module ecore = m.def_submodule("core", "Core errors");
    PyErrorsCore(ecore);
}
