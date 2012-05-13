
include /usr/local/trilinos/include/Makefile.export.Anasazi
include /usr/local/trilinos/include/Makefile.export.AztecOO
include /usr/local/trilinos/include/Makefile.export.Ifpack

CXX = $(Anasazi_CXX_COMPILER)
CXXFLAGS = -O0 -g -Wall

INCLUDE = -I/usr/local/trilinos/include -I/ports/include -I/data/software/include -I.
LINK = -L/usr/local/trilinos/lib -L/data/software/lib $(Anasazi_LIBRARIES) $(AztecOO_LIBRARIES) $(Ifpack_LIBRARIES) -lepetraext -framework vecLib



all : test_anasazi test_anasazi_shiftinvert

test_anasazi_shiftinvert : test_anasazi_shiftinvert.o test_anasazi_tools.o
	$(CXX) -o test_anasazi_shiftinvert test_anasazi_shiftinvert.o test_anasazi_tools.o $(LINK)

test_anasazi : test_anasazi.o test_anasazi_tools.o
	$(CXX) -o test_anasazi test_anasazi.o test_anasazi_tools.o $(LINK)

%.o : %.cpp test_anasazi.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

install : test_anasazi test_anasazi_shiftinvert
	cp test_anasazi test_anasazi_shiftinvert /data/software/bin/

clean :
	rm -f *.o test_anasazi test_anasazi_shiftinvert *~

