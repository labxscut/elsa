#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "compcore.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_compcore_gpu, m) {
    m.doc() = "GPU backend module for LSA/LLA compatibility; it shares the same interface as lsa._compcore.";

    m.def("test", &test);
    m.def("calc_LA", &calc_LA,
          py::arg("x"),
          py::arg("y"),
          py::arg("z"));

    py::class_<LSA_Data>(m, "LSA_Data")
        .def(py::init<>())
        .def(py::init<int, VectorDouble, VectorDouble>())
        .def_readwrite("max_shift", &LSA_Data::max_shift)
        .def_readwrite("X", &LSA_Data::X)
        .def_readwrite("Y", &LSA_Data::Y)
        .def("assign", &LSA_Data::assign);

    py::class_<LSA_Result>(m, "LSA_Result")
        .def(py::init<>())
        .def_readwrite("score", &LSA_Result::score)
        .def_readwrite("trace", &LSA_Result::trace);

    m.def("DP_lsa", &DP_lsa,
          py::arg("data"),
          py::arg("keep_trace") = true);

    py::class_<LLA_Data>(m, "LLA_Data")
        .def(py::init<>())
        .def(py::init<int, VectorDouble, VectorDouble, VectorDouble>())
        .def_readwrite("max_shift", &LLA_Data::max_shift)
        .def_readwrite("X", &LLA_Data::X)
        .def_readwrite("Y", &LLA_Data::Y)
        .def_readwrite("Z", &LLA_Data::Z)
        .def("assign", &LLA_Data::assign);

    py::class_<LLA_Result>(m, "LLA_Result")
        .def(py::init<>())
        .def_readwrite("score", &LLA_Result::score)
        .def_readwrite("trace", &LLA_Result::trace);

    m.def("DP_lla", &DP_lla,
          py::arg("data"),
          py::arg("keep_trace") = true);

    m.attr("__version__") = "2.0.0";
    m.attr("__author__") = "Li Charles Xia";
    m.attr("__copyright__") = "Copyright (c) 2008-2024 Li Charles Xia";
    m.attr("__license__") = "BSD";
}