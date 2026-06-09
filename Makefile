# Heat Equation — build the serial and MPI solvers into bin/
#
#   make            # build both binaries
#   make serial     # build only bin/heat-big
#   make mpi        # build only bin/heat-mpi-big
#   make clean      # remove bin/

CC      ?= cc
MPICC   ?= mpicc
CFLAGS  ?= -O3 -Wall -Wextra
LDLIBS  ?= -lm

SRC := src
BIN := bin

SERIAL := $(BIN)/heat-big
MPI    := $(BIN)/heat-mpi-big

.PHONY: all serial mpi clean
all: serial mpi

serial: $(SERIAL)
mpi:    $(MPI)

$(SERIAL): $(SRC)/heat-big.c | $(BIN)
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

$(MPI): $(SRC)/heat-mpi-big.c | $(BIN)
	$(MPICC) $(CFLAGS) $< -o $@ $(LDLIBS)

$(BIN):
	mkdir -p $(BIN)

clean:
	rm -rf $(BIN)
