# export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5
# export ASAN_OPTIONS=symbolize=1

CXX = clang++
CXXFLAGS = -W -Wall -Wextra -std=c++11 -pthread -O3 -march=native -mtune=native
OBJECTS = *.o
CIPHER = ../minikiasu/tablebased

ifdef DEBUG
CXXFLAGS += -g3 -DDEBUG -fsanitize=undefined -fsanitize=address -fsanitize=alignment -ftrapv -fno-omit-frame-pointer -fno-optimize-sibling-calls
endif

.PHONY: all clean format test-kiasu boomerang-differentials create-kiasuh gen-tables 

all: boomerang-differentials

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(CIPHER) -c $< -o $@

$(CIPHER)/%.o: $(CIPHER)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(CIPHER) -c $< -o $@

boomerang-differentials: create-kiasuh main.o $(CIPHER)/kiasu.o
	$(CXX) $(CXXFLAGS) -I$(CIPHER) main.cpp $(CIPHER)/kiasu.cpp -o $@

test-kiasu: create-kiasuh $(CIPHER)/kiasu.o
	$(CXX) $(CXXFLAGS) -I$(CIPHER) $(CIPHER)/test-kiasu7.cpp $(CIPHER)/kiasu.cpp -o $@

create-kiasuh:
	cp -f $(CIPHER)/kiasu7.h $(CIPHER)/kiasu.h

format:
	astyle *.cpp *.h

clean:
	rm -f *.o bench core architectures boomerang-differentials test-kiasu

