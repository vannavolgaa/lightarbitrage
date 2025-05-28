#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_errors_frameworks.h"
#include "py_errors_math.h"
#include "../../../../src/pybindings/errors/core/datastructure/py_errors_datastructure.h"
#include "py_errors_assetpricer.h"
#include "py_errors_calibration.h"

namespace py = pybind11;

void PyErrorsCore(py::module_& m)
{
    py::module emath = m.def_submodule("math", "Math errors");
    py::module eframeworks = m.def_submodule("frameworks", "Frameworks errors");
    py::module edts = m.def_submodule("datastructures", "Datastructure errors");
    py::module eassetpricer = m.def_submodule("assetpricer", "Asset pricer errors");
    py::module ecalibration = m.def_submodule("calibration", "Calibration errors");
    
    PyErrorsFrameworks(eframeworks);
    PyErrorsMath(emath);
    PyErrorsDataStructure(edts);
    PyErrorsAssetPricer(eassetpricer);
    PyErrorsCalibration(ecalibration);
    
}