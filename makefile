CXXFLAGS = -O3 -std=c++11
#CXXFLAGS = -g -std=c++11 # for debug
BIN_DIR = bin/

all:
	mkdir -p $(BIN_DIR)
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

PREFIX = /usr/local
EXEC_PREFIX= $(PREFIX)
EXEC_BIN_DIR = $(EXEC_PREFIX)/bin/

install: all
	mkdir -p $(EXEC_BIN_DIR)
	cp $(BIN_DIR)* $(EXEC_BIN_DIR)

clean:
	@cd src && $(MAKE) clean