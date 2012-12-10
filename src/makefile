# Generic makefile for flowVC
#
#
ifeq ($(mode),debug)
   CFLAGS = -g -Wall -O0   
else
   mode = release
   CFLAGS = -Wall -O3 
endif
LFLAGS = -lm
CC = gcc
SRC = exposuretime.c ftle.c globals.c integration.c io.c main.c memory.c mesh.c mymath.c parameters.c \
	residencetime.c strainrate.c tracers.c velocity.c velout.c vorticity.c
OBJ = $(SRC:.c=.o)
EXE = flowVC

.PHONY: all
all: info $(EXE)

info:
ifneq ($(mode),release)
ifneq ($(mode),debug)
    @echo "Invalid build mode." 
    @echo "Please use 'make mode=release' or 'make mode=debug'"
    @exit 1
endif
endif
	@echo "Building in "$(mode)" mode..."

$(EXE): $(OBJ)
	$(CC) $(LFLAGS) $(OBJ) -o $@

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -vf *.o $(EXE)


