# Makefile for building an example plugin

CXX = $(shell harpconfig --cxx)
CXXFLAGS = $(shell harpconfig --cxxflags --cppflags) -I.
PLUG_FLAGS = $(shell harpconfig --plugflags)
PLUG_LINK = $(shell harpconfig --pluglink)
LINK = $(shell harpconfig --link)

# the output plugin files *MUST* have names of the format:
#
# harp_plugin_(type)_(name).so
#
# where type = {image, spec, psf}

all : harp_plugin_example.so harp_example

harp_plugin_example.so : harp_plugin_example.o
	$(CXX) $(PLUG_FLAGS) -o $@ $< $(PLUG_LINK) $(LINK)

harp_plugin_example.o : harp_plugin_example.cpp harp_example.hpp
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -o $@ -c $<

harp_example : harp_example.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINK)

harp_example.o : harp_example.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<


clean :
	rm -f *.so *.o harp_example *~

