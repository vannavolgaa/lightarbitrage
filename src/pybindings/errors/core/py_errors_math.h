#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/errors/core/math.h"

namespace py = pybind11;

void PyErrorsMath(pybind11::module_& m)
{
    using Code = MathErrorCode;
    using Error = MathError;

    py::enum_<Code>(m, "MathErrorCode")
        .value("NEWTON_RAPHSON_DERIVATIVE_CLOSE_TO_ZERO", MathErrorCode::NEWTON_RAPHSON_DERIVATIVE_CLOSE_TO_ZERO)
        .value("INTERPOLATION_OUT_OF_RANGE", MathErrorCode::INTERPOLATION_OUT_OF_RANGE)
        .value("INTERPOLATION_WRONG_VECTOR_SIZE", MathErrorCode::INTERPOLATION_WRONG_VECTOR_SIZE)
        .value("NORMAL_DISTRIBUTION_INVALID_SIGMA", MathErrorCode::NORMAL_DISTRIBUTION_INVALID_SIGMA)
        .value("MISATCH_VECTOR_SIZE_LOSS_FUNCTION", MathErrorCode::MISATCH_VECTOR_SIZE_LOSS_FUNCTION)
        .value("LOSS_ZERO_VALUE_IN_TRUE_VALUES", MathErrorCode::LOSS_ZERO_VALUE_IN_TRUE_VALUES)
        .value("LOSS_NEGATIVE_VALUE", MathErrorCode::LOSS_NEGATIVE_VALUE)
        .value("LOSS_UNKNOWN_TYPE", MathErrorCode::LOSS_UNKNOWN_TYPE)
        .value("REGRESSION_INVALID_X_MATRIX", MathErrorCode::REGRESSION_INVALID_X_MATRIX)
        .value("REGRESSION_MISMATCH_X_Y_SIZE", MathErrorCode::REGRESSION_MISMATCH_X_Y_SIZE)
        .value("UNIFORM_INVALID_BOUNDS", MathErrorCode::UNIFORM_INVALID_BOUNDS)
        .value("CORRELATION_MATRIX_NOT_SQUARE", MathErrorCode::CORRELATION_MATRIX_NOT_SQUARE)
        .value("CORRELATION_MATRIX_NOT_SYMMETRIC", MathErrorCode::CORRELATION_MATRIX_NOT_SYMMETRIC)
        .value("CORRELATION_MATRIX_NOT_POSITIVE_DEFINITE", MathErrorCode::CORRELATION_MATRIX_NOT_POSITIVE_DEFINITE)
        .value("CORRELATION_INVALID_COEFFICIENT", MathErrorCode::CORRELATION_INVALID_COEFFICIENT)
        .value("CORRELATION_INVALID_COEFFICIENT_SIZE", MathErrorCode::CORRELATION_INVALID_COEFFICIENT_SIZE)
        .export_values();

    // Exception translation
    static py::exception<MathError> py_error(m, "MathError");

    py::class_<MathError>(m, "ArbitrageMathError")
        .def(py::init<MathErrorCode>(), py::arg("code"))
        .def("get_code", &MathError::get_code)
        .def("__str__", [](const MathError& e) { return e.what(); });

    // Register exception translator
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const MathError& e) {
            PyErr_SetString(py_error.ptr(), e.what());
        }
    });

}



