# Makefile for building an example plugin

MPICXX = $(shell harpconfig --mpicxx)

CXX = $(shell harpconfig --cxx)
CXXFLAGS = $(shell harpconfig --cxxflags --cppflags) -I.

LINK = $(shell harpconfig --link)
MPILINK = $(shell harpconfig --mpilink)

PLUG_FLAGS = $(shell harpconfig --plugflags)
PLUG_LINK = $(shell harpconfig --pluglink)
PLUG_EXT = $(shell harpconfig --plugext)

# the output plugin files *MUST* have names of the format:
#
# harp_plugin_(type)_(name).$(PLUG_EXT)
#
# where type = {image, spec, psf}

all : harp_example harp_mpi_example


harp_plugin_example.$(PLUG_EXT) : harp_plugin_example.o
	$(CXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< $(LINK)

harp_plugin_example.o : harp_plugin_example.cpp harp_example.hpp
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -o $@ -c $<


harp_example : harp_example.o harp_plugin_example.$(PLUG_EXT)
	$(CXX) -o $@ $< $(LINK)

harp_example.o : harp_example.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<


harp_mpi_example : harp_mpi_example.o harp_plugin_example.$(PLUG_EXT)
	$(MPICXX) -o $@ $< $(MPILINK)

harp_mpi_example.o : harp_mpi_example.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ -c $<


clean :
	rm -f *.$(PLUG_EXT) *.o harp_example harp_mpi_example *~

