#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../src/pybindings/core/py_core.h"
#include "../../src/pybindings/errors/py_errors_arbitrage.h"

namespace py = pybind11;

PYBIND11_MODULE(lightarbitrage, m) {
    
    m.doc() = "This is light Arbitrage, written in C++ and wrapped for python";
    
    py::module core = m.def_submodule("core", "Definition of Arbitrage's core module.");
    py::module errors = m.def_submodule("errors", "Definition of Arbitrage's error handling module.");
    PyErrorsArbitrage(errors);
    PyCore(core);
    
}