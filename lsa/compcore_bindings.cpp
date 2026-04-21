#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "compcore.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_compcore, m) {
    m.doc() = "Python bindings for LSA (Local Similarity Analysis) and LA/LLA (Liquid Association) algorithms";

    // ====================================
    // Static LA (Liquid Association) Functions
    // ====================================
    
    // Bind basic test function
    m.def("test", &test, 
          "Simple test function that prints 'tested'");

    // Bind static LA calculation function
    m.def("calc_LA", &calc_LA,
          py::arg("x"),  // First time series
          py::arg("y"),  // Second time series
          py::arg("z"),  // Third time series (mediator)
          "Calculate static Liquid Association score.\n"
          "Args:\n"
          "    x: First time series data\n"
          "    y: Second time series data\n"
          "    z: Third time series data (mediator)\n"
          "Returns:\n"
          "    float: LA score between the three series");

    // ====================================
    // LSA (Local Similarity Analysis) Bindings
    // ====================================
    
    // Bind LSA_Data class - holds input data for LSA computation
    py::class_<LSA_Data>(m, "LSA_Data")
        // Default constructor
        .def(py::init<>())
        // Constructor with parameters (delay limit and two time series)
        .def(py::init<int, VectorDouble, VectorDouble>(),
             py::arg("max_shift"),   // Maximum allowed delay between series
             py::arg("X"),           // First time series
             py::arg("Y"))           // Second time series
        // Expose member variables for direct access from Python
        .def_readwrite("max_shift", &LSA_Data::max_shift, "Maximum allowed delay between series")
        .def_readwrite("X", &LSA_Data::X, "First time series data")
        .def_readwrite("Y", &LSA_Data::Y, "Second time series data")
        .def("assign", &LSA_Data::assign,
             py::arg("max_shift"),
             py::arg("X"),
             py::arg("Y"),
             R"pbdoc(
             Assign new values to the LSA_Data object.

             Args:
                 max_shift (int): Maximum allowed delay between series
                 X (list): First time series data
                 Y (list): Second time series data

             Example:
                 >>> lsa_data = LSA_Data()
                 >>> lsa_data.assign(2, [1.0, 2.0, 3.0], [2.0, 3.0, 4.0])
             )pbdoc"
        );

    // Bind LSA_Result class - holds results of LSA computation
    py::class_<LSA_Result>(m, "LSA_Result")
        .def(py::init<>())
        .def_readwrite("score", &LSA_Result::score, "LSA similarity score")
        .def_readwrite("trace", &LSA_Result::trace, "Alignment trace for best score");

    // Bind DP_lsa function - performs LSA computation
    m.def("DP_lsa", &DP_lsa, 
          py::arg("data"),           // Input LSA_Data object
          py::arg("keep_trace") = true,  // Whether to compute alignment trace
          "Compute Local Similarity Analysis (LSA) score using dynamic programming.\n"
          "Args:\n"
          "    data: LSA_Data object containing input series\n"
          "    keep_trace: If True, compute alignment trace\n"
          "Returns:\n"
          "    LSA_Result object containing score and optional trace");

    // ====================================
    // LLA (Local Liquid Association) Bindings
    // ====================================
    
    // Bind LLA_Data class - holds input data for LLA computation
    py::class_<LLA_Data>(m, "LLA_Data")
        // Default constructor
        .def(py::init<>())
        // Constructor with parameters (delay limit and three time series)
        .def(py::init<int, VectorDouble, VectorDouble, VectorDouble>(),
             py::arg("max_shift"),   // Maximum allowed delay
             py::arg("X"),           // First time series
             py::arg("Y"),           // Second time series
             py::arg("Z"))           // Third time series (mediator)
        // Expose member variables for direct access from Python
        .def_readwrite("max_shift", &LLA_Data::max_shift, "Maximum allowed delay between series")
        .def_readwrite("X", &LLA_Data::X, "First time series data")
        .def_readwrite("Y", &LLA_Data::Y, "Second time series data")
        .def_readwrite("Z", &LLA_Data::Z, "Third time series data (mediator)")
        .def("assign", &LLA_Data::assign,
             py::arg("max_shift"),
             py::arg("X"),
             py::arg("Y"),
             py::arg("Z"),
             R"pbdoc(
             Assign new values to the LLA_Data object.

             Args:
                 max_shift (int): Maximum allowed delay between series
                 X (list): First time series data
                 Y (list): Second time series data
                 Z (list): Third time series data (mediator)

             Example:
                 >>> lla_data = LLA_Data()
                 >>> lla_data.assign(2, [1.0, 2.0], [2.0, 3.0], [3.0, 4.0])
             )pbdoc"
        );

    // Bind LLA_Result class - holds results of LLA computation
    py::class_<LLA_Result>(m, "LLA_Result")
        .def(py::init<>())
        .def_readwrite("score", &LLA_Result::score, "LLA score")
        .def_readwrite("trace", &LLA_Result::trace, "Alignment trace for best score");

    // Bind DP_lla function - performs LLA computation
    m.def("DP_lla", &DP_lla,
          py::arg("data"),           // Input LLA_Data object
          py::arg("keep_trace") = true,  // Whether to compute alignment trace
          "Compute Local Liquid Association (LLA) score using dynamic programming.\n"
          "Args:\n"
          "    data: LLA_Data object containing input series\n"
          "    keep_trace: If True, compute alignment trace\n"
          "Returns:\n"
          "    LLA_Result object containing score and optional trace");

    // ====================================
    // Module Metadata
    // ====================================
    m.attr("__version__") = "2.0.0";
    m.attr("__author__") = "Li Charles Xia";
    m.attr("__copyright__") = "Copyright (c) 2008-2024 Li Charles Xia";
    m.attr("__license__") = "BSD";
} 