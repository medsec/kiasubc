# export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5
# export ASAN_OPTIONS=symbolize=1

CXX = clang++
CXXFLAGS = -W -Wall -Wextra -std=c++11 -pthread -march=native -maes -mavx -O3
OBJECTS = *.o
CIPHER = ../kiasu

ifdef DEBUG
CXXFLAGS += -g3 -DDEBUG -fsanitize=undefined -fsanitize=address -fsanitize=alignment -ftrapv -fno-omit-frame-pointer -fno-optimize-sibling-calls
endif

.PHONY: all clean format test-kiasu boomerang-differentials

all: boomerang-differentials

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(CIPHER) -c $< -o $@

$(CIPHER)/%.o: $(CIPHER)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(CIPHER) -c $< -o $@

boomerang-differentials: main.o $(CIPHER)/kiasu.o
	$(CXX) $(CXXFLAGS) -I$(CIPHER) main.cpp $(CIPHER)/kiasu.cpp -o $@

test-kiasu: $(CIPHER)/kiasu.o
	$(CXX) $(CXXFLAGS) -I$(CIPHER) $(CIPHER)/test-kiasu.cpp $(CIPHER)/kiasu.cpp -o $@

format:
	astyle *.cpp *.h

clean:
	rm -f *.o bench core architectures boomerang-differentials test-kiasu

