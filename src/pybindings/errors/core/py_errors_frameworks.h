#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/errors/core/frameworks.h"

namespace py = pybind11;

void PyErrorsFrameworks(py::module_& m)
{
    using Code = FrameworksErrorCode;
    using Error = FrameworksError;

    // Bind enum
    py::enum_<Code>(m, "FrameworksErrorCode")
        .value("INVALID_NORMALIZED_BLACK_PRICE", Code::INVALID_NORMALIZED_BLACK_PRICE)
        .value("JUMP_WINGS_SVI_INVALID_PARAM", Code::JUMP_WINGS_SVI_INVALID_PARAM)
        .value("RAW_SVI_INVALID_PARAM_B", Code::RAW_SVI_INVALID_PARAM_B)
        .value("RAW_SVI_INVALID_PARAM_P", Code::RAW_SVI_INVALID_PARAM_P)
        .value("RAW_SVI_INVALID_PARAM_BETA", Code::RAW_SVI_INVALID_PARAM_BETA)
        .value("RAW_SVI_INVALID_PARAM_S", Code::RAW_SVI_INVALID_PARAM_S)
        .value("RAW_SVI_CONDITION", Code::RAW_SVI_CONDITION)
        .value("SSVI_INVALID_PARAM", Code::SSVI_INVALID_PARAM)
        .value("VASICEK_INVALID_PARAM", Code::VASICEK_INVALID_PARAM)
        .value("BLACK_SCHOLES_MODEL_INVALID_PARAM", Code::BLACK_SCHOLES_MODEL_INVALID_PARAM)
        .value("BARONE_ADESI_WHALEY_MODEL_INVALID_PARAM", Code::BARONE_ADESI_WHALEY_MODEL_INVALID_PARAM)
        .value("LETS_BE_RATIONAL_INVALID_BS_PRICE", Code::LETS_BE_RATIONAL_INVALID_BS_PRICE)
        .value("HESTON_INVALID_PARAM", Code::HESTON_INVALID_PARAM)
        .value("HESTON_UNKNOWN_DISCRETIZATION_METHOD", Code::HESTON_UNKNOWN_DISCRETIZATION_METHOD)
        .export_values();
    
    static py::exception<Error> py_error(m, "FrameworksError");

    py::class_<Error>(m, "ArbitrageFrameworksError")
        .def(py::init<Code>(), py::arg("code"))
        .def("get_code", &Error::get_code)
        .def("__str__", [](const Error& e) { return e.what(); });

    // Register exception translator
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const Error& e) {
            PyErr_SetString(py_error.ptr(), e.what());
        }
    });


};
