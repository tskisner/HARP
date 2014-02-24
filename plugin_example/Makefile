# Makefile for building an example plugin

MPICXX = $(shell harpconfig --mpicxx)
CXX = $(shell harpconfig --cxx)
CXXFLAGS = $(shell harpconfig --cxxflags --cppflags) -I.
PLUG_FLAGS = $(shell harpconfig --plugflags)
PLUG_LINK = $(shell harpconfig --pluglink)
LINK = $(shell harpconfig --link)
MPILINK = $(shell harpconfig --mpilink)

# the output plugin files *MUST* have names of the format:
#
# harp_plugin_(type)_(name).so
#
# where type = {image, spec, psf}

all : harp_example harp_mpi_example


harp_plugin_example.so : harp_plugin_example.o
	$(CXX) $(PLUG_FLAGS) -o $@ $< $(PLUG_LINK) $(LINK)

harp_plugin_example.o : harp_plugin_example.cpp harp_example.hpp
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -o $@ -c $<


harp_mpi_plugin_example.so : harp_mpi_plugin_example.o
	$(MPICXX) $(PLUG_FLAGS) -o $@ $< $(PLUG_LINK) $(MPILINK)

harp_mpi_plugin_example.o : harp_plugin_example.cpp harp_example.hpp
	$(MPICXX) $(CXXFLAGS) $(PLUG_FLAGS) -DUSE_MPI -o $@ -c $<


harp_example : harp_example.o harp_plugin_example.so
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINK)

harp_example.o : harp_example.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<


harp_mpi_example : harp_mpi_example.o harp_mpi_plugin_example.so
	$(MPICXX) $(CXXFLAGS) -o $@ $< $(MPILINK)

harp_mpi_example.o : harp_mpi_example.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ -c $<


clean :
	rm -f *.so *.o harp_example harp_mpi_example *~

