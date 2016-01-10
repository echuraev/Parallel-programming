SOURCES=main.cpp
TARGET=swe
CXX=mpiCC
CXX_FLAGS=

.PHONY: $(TARGET)

all: $(TARGET)

$(TARGET): $(SOURCES) 
	$(CXX) $(CXX_FLAGS) $^ -o $@

.PHONY: clean

clean: 
	rm -vf $(TARGET) 
