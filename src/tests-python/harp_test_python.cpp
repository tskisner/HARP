
#include <iostream>

#include <harp.hpp>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

using namespace std;
using namespace harp;

namespace py = boost::python


int add_five(int x) {
  return x + 5;
}

BOOST_PYTHON_MODULE(Pointless)
{
  py::def("add_five", add_five);
}


int main ( int argc, char *argv[] ) {

  Py_Initialize();

  try {

    initPointless();

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import harp_test_python_module", main_namespace );
    py::object harptest = py::eval ( "harp_test_python_module.HarpTest(\"my/path\", 1.5, 4.0)", main_namespace );

    string path = py::extract < string > ( harptest.attr("get_path")() );

    cout << "path = " << path << endl;

    double diff = py::extract < double > ( harptest.attr("get_diff")() );

    cout << "diff = " << diff << endl;

    py::exec ( "import Pointless", main_namespace );
    py::object pntless = py::eval ( "Pointless.add_five(4)", main_namespace );
    int result = py::extract < int > ( pntless );

    cout << "pointless = " << result << endl;

  } catch ( py::error_already_set ) {
    PyErr_Print();
  }

  // boost python docs recommend NOT calling this...
  //Py_Finalize();

  return 0;
}
