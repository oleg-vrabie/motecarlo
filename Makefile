SOURCES=$(HOME)/devel/sources
vpath %.cc $(SOURCES) 
vpath %.cpp $(SOURCES)
vpath %.c $(SOURCES)

CXX = g++
CXXFLAGS = -fopenmp -O3 -std=c++11 -I. -I $(SOURCES)/lib
LIBS = -lgsl -lgslcblas 

PROGRAM = montecarlo
OBJS = montecarlo.o ex1.o

$(PROGRAM): $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LIBS)
	
clean:
	rm -f $(PROGRAMS) $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 
	
%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 
	
%.o: %.c $(HDRS)
	gcc $(CXXFLAGS) -c $< 
