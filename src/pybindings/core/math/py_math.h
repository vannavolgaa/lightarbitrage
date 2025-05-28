#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "py_interpolation2D.h"
#include "py_optimization.h"
#include "py_quadratures.h"
#include "py_erfcody.h"
#include "py_lossfunction.h"
#include "py_regression.h"
#include "../../../../src/pybindings/core/math/probability/py_probabillity.h"

void PyMath(py::module_& m)
{
    py::module linearinterpolation2D = m.def_submodule("linearinterpolation2D", "Definition of 2D linear interpolation objects");
    py::module quadratures = m.def_submodule("quadratures", "Definition of quadratures");
    py::module optimization = m.def_submodule("optimization", "Definition of optimization objects");
    py::module probability = m.def_submodule("probability", "Definition of probability objects");
    py::module erfcody = m.def_submodule("erfcody", "Definition of erfcody objects");
    py::module lossfunctions = m.def_submodule("lossfunctions", "Definition of loss functions");
    py::module regression = m.def_submodule("regression", "Definition of regression objects");

    PyErfCody(erfcody);
    PyProbability(probability);
    PyInterpolation2D(linearinterpolation2D);
    PyQuadratures(quadratures);
    PyOptimization(optimization);
    PyLossFunctions(lossfunctions);
    PyRegression(regression);

}