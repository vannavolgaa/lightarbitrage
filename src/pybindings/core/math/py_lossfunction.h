#pragma once 
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/core/math/lossfunction.h"

namespace py = pybind11;

void PyLossFunctions(py::module_& m)
{
    py::enum_<LossType>(m, "LossType")
        .value("MSE", LossType::MSE)
        .value("RMSE", LossType::RMSE)
        .value("MAE", LossType::MAE)
        .value("MAPE", LossType::MAPE)
        .value("MSLE", LossType::MSLE)
        .value("MSPE", LossType::MSPE)
        .value("RMPSE", LossType::RMPSE)
        .export_values();

    // Bind the compute_loss function
    m.def("compute_loss", &compute_loss, py::arg("estimate"), py::arg("true_values"), py::arg("loss_type"),
          "Compute the loss between estimated and true values using the specified loss type.");
}