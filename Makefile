# Compiler and flags
CC = gcc
CFLAGS = -Ofast -lm

# Source files and output executable
SRC = main.c PSO.c OF.c
OBJ = $(SRC:.c=.o)
EXEC = pso

# Default target
all: $(EXEC)

# Build the executable
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(CFLAGS)

# Compile source files into object files
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Clean up object files and the executable
clean:
	rm -f $(OBJ) $(EXEC)
