SOURCES=sources/sharp.cpp sources/function.cpp sources/parse_conf.cpp 
BIN_DIR=bin/
INPUT_DIR=input/
TMPL_DIR=sources/tmpl/
INPUT_CONF=$(INPUT_DIR)input.conf
TARGET=sharp
CXX=g++
CXX_FLAGS=-g

.PHONY: $(TARGET)

all: project_conf $(TARGET)

$(TARGET): $(SOURCES) 
	$(CXX) $(CXX_FLAGS) $^ -o $(BIN_DIR)$@

project_conf: $(BIN_DIR) $(INPUT_DIR) $(INPUT_CONF) 

$(BIN_DIR) $(INPUT_DIR):
	@if [ ! -d $@ ]; then mkdir $@; fi

$(INPUT_CONF): $(TMPL_DIR)input.conf
	@cp -f $(TMPL_DIR)* $(INPUT_DIR)

.PHONY: clean

clean: 
	rm -vf $(BIN_DIR)$(TARGET) 
