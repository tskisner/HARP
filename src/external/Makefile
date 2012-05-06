
include /usr/local/trilinos/include/Makefile.export.Anasazi

CXX = $(Anasazi_CXX_COMPILER)
CXXFLAGS = -O0 -g -Wall

INCLUDE = -I/usr/local/trilinos/include -I/ports/include
LINK = -L/usr/local/trilinos/lib $(Anasazi_LIBRARIES) -lepetraext -framework vecLib



all : test_anasazi

test_anasazi : test_anasazi.o
	$(CXX) -o test_anasazi test_anasazi.o $(LINK)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

install : test_anasazi
	cp test_anasazi /data/software/bin/

clean :
	rm -f *.o test_anasazi *~

