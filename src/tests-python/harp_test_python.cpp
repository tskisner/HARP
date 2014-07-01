
#include <iostream>

#include <harp.hpp>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

using namespace std;
using namespace harp;


int add_five(int x) {
  return x + 5;
}


BOOST_PYTHON_MODULE(Pointless)
{
  boost::python::def("add_five", add_five);
}


int main ( int argc, char *argv[] ) {

  Py_Initialize();

  try {

    initPointless();

    boost::python::object main_module = boost::python::import("__main__");
    boost::python::object main_namespace = main_module.attr("__dict__");

    boost::python::object ignored = boost::python::exec(
      "import Pointless\n"
      "print Pointless.add_five(4)\n"
      "result = 5 ** 2\n"
      "print result\n"
      , main_namespace);

  } catch ( boost::python::error_already_set ) {
    PyErr_Print();
  }

  // boost python docs recommend NOT calling this...
  //Py_Finalize();

  return 0;
}
