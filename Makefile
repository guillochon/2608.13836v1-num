CXX ?= g++
OMP ?= -fopenmp
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -march=native
LDFLAGS ?=

TARGET := brenti_search
SRC := src/main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(OMP) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
