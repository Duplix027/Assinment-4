#ifndef UTILITY_H
#define UTILITY_H

// Function pointer type for objective functions
typedef double (*ObjectiveFunction)(int, double *);

// Structure for bounds of each variable
typedef struct Bound {
    double lowerBound; // The lowr limit for a variable in the search space
    double upperBound; // The upper limit for a variable in the search space
} Bound;

// Structure for a single particle
typedef struct Particle {
    double *position;       // Current postion of the particle in the serch space
    double *velocity;       // Current velocioty or direction of movement
    double *best_position;  // The personal best postion of this particle
    double value;           // Fitness value of the current postion
    double best_value;      // Best fitness value ever achieved by the particle
} Particle;

// Function prototypes

// Helper function to generate a random double between min and max.
// NOTE: Random values are used to create diversitiy in particle movement
double random_double(double min, double max);

// Function to initialize all particles
// Sets their positions within the search space bounds and assigns small random velocities
// The glboal best position and value are also updated if the initialized particle has better fitness

//CHATGPT was used to help come up with function parameters
void create_p(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds, ObjectiveFunction objective_function, double *global_best_position, double *global_best_value);

// Updates the velocities and positions of particles
// Adjusts based on personal and global bests using the PSO formula
// Also clamps the velocities to prevent "too fast" movements, which could break convergence
// CHATGPT was used to help come up with this idea of clamping velocities

//CHATGPT was used to help come up with function parameters
void repos_p(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds,
                      ObjectiveFunction objective_function, double *global_best_position, double *global_best_value, double w, double c1, double c2);

// Free the dynamically allocated memory for particle attributes
void free_particles(Particle *particles, int NUM_PARTICLES);

// The main PSO function
// Orchestrates the whole process:
// - Initializes particles
// - Iteratively updates particles' velocities and positions
// - Tracks the global best solution
// Returns the fitness of the best solution found after all iterations
double pso(ObjectiveFunction objective_function, int NUM_VARIABLES, Bound *bounds, int NUM_PARTICLES, int MAX_ITERATIONS, double *best_position);


#endif 
