
include /usr/local/trilinos/include/Makefile.export.Anasazi
include /usr/local/trilinos/include/Makefile.export.AztecOO
include /usr/local/trilinos/include/Makefile.export.Belos
include /usr/local/trilinos/include/Makefile.export.Ifpack

CXX = $(Anasazi_CXX_COMPILER)
CXXFLAGS = -O3 -Wall -fopenmp
#CXXFLAGS = -O0 -g -Wall -fopenmp

INCLUDE = -I/usr/local/trilinos/include -I/ports/include -I/data/software/include -I.
LINK = -L/usr/local/trilinos/lib -L/data/software/lib $(Anasazi_LIBRARIES) $(AztecOO_LIBRARIES) $(Belos_LIBRARIES) $(Ifpack_LIBRARIES) -lepetraext -framework vecLib -fopenmp



all : test_anasazi_shiftinvert test_anasazi test_anasazi_vis

test_anasazi_shiftinvert : test_anasazi_shiftinvert.o test_anasazi_tools.o
	$(CXX) -o test_anasazi_shiftinvert test_anasazi_shiftinvert.o test_anasazi_tools.o $(LINK)

test_anasazi : test_anasazi.o test_anasazi_tools.o
	$(CXX) -o test_anasazi test_anasazi.o test_anasazi_tools.o $(LINK)

test_anasazi_vis : test_anasazi_vis.o test_anasazi_tools.o
	$(CXX) -o test_anasazi_vis test_anasazi_vis.o test_anasazi_tools.o $(LINK)

%.o : %.cpp test_anasazi.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

install : test_anasazi test_anasazi_shiftinvert test_anasazi_vis
	cp test_anasazi test_anasazi_shiftinvert test_anasazi_vis /data/software/bin/

clean :
	rm -f *.o test_anasazi test_anasazi_shiftinvert test_anasazi_vis *~

