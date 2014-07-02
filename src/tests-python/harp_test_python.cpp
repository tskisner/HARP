
#include <iostream>

#include <harp.hpp>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

using namespace std;
using namespace harp;

namespace py = boost::python;


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

    string cpath = "/my/favorite/path";
    double minval = 1.5;
    double maxval = 4.0;

    ostringstream com;
    
    com.str("");
    com << "test = harp_test_python_module.HarpTest(\"" << cpath << "\", " << minval << ", " << maxval << ")";

    // instantiate test in the main dictionary
    py::exec ( com.str().c_str(), main_namespace );

    // get a handle to test in c++
    py::object harptest = py::eval ( "test", main_namespace );

    // pickle this and return a string to c++
    py::exec ( "import pickle", main_namespace );
    string pickle = py::extract < string > ( py::eval ( "pickle.dumps( test, -1 )", main_namespace ) );

    // now unpickle to a new instance in python
    main_module.attr ( "picklestr" ) = pickle;
    py::exec ( "test2 = pickle.loads( picklestr )", main_namespace );

    py::object harptest2 = py::eval ( "test2", main_namespace );

    string path = py::extract < string > ( harptest.attr("get_path")() );
    cout << "1 : path = " << path << endl;

    double diff = py::extract < double > ( harptest.attr("get_diff")() );
    cout << "1 : diff = " << diff << endl;

    path = py::extract < string > ( harptest2.attr("get_path")() );
    cout << "2 : path = " << path << endl;

    diff = py::extract < double > ( harptest2.attr("get_diff")() );
    cout << "2 : diff = " << diff << endl;

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
