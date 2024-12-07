# Compiler and flags
CC = gcc
CFLAGS = -Ofast -lm

# Source files and output executable
SRC = main.c PSO.c OF.c
EXEC = pso

# Default target
all: $(EXEC)

# Build the executable directly from source files
$(EXEC):
	$(CC) $(SRC) -o $(EXEC) $(CFLAGS)

# Clean up the executable
clean:
	rm -f $(EXEC)
