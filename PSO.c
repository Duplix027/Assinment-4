#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "utility.h"

// Helper function to generate random numbers in a range
// This function generates a double precision number between `min` and `max`.
// It is used to randomize particle positions and velocities.
double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Function to initialize particles
// This function sets up the initial state of all particles.
// Positions are randomly assigned within the provided bounds.
// Velocities are also randomized but kept small to prevent erratic behavior.
// The global best position is updated during initialization if a particle achieves a better fitness.
void initialize_particles(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds, ObjectiveFunction objective_function, double *global_best_position, double *global_best_value) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].position = malloc(NUM_VARIABLES * sizeof(double)); // Allocate memory for position
        particles[i].velocity = malloc(NUM_VARIABLES * sizeof(double)); // Allocate memory for velocity
        particles[i].best_position = malloc(NUM_VARIABLES * sizeof(double)); // Allocate memory for best position
        particles[i].best_value = DBL_MAX; // Set the best value to a very large number initially --> CHATGPT suggested this

        for (int j = 0; j < NUM_VARIABLES; j++) {
            particles[i].position[j] = random_double(bounds[j].lowerBound, bounds[j].upperBound); // Random position
            particles[i].velocity[j] = random_double(-1.0, 1.0); // Random velocity in a small range
            particles[i].best_position[j] = particles[i].position[j]; // Initialize personal best position
        }

        particles[i].value = objective_function(NUM_VARIABLES, particles[i].position); // Evaluate the fitness
        if (particles[i].value < *global_best_value) { // Update global best if this particle is better
            *global_best_value = particles[i].value;
            for (int j = 0; j < NUM_VARIABLES; j++) {
                global_best_position[j] = particles[i].position[j];
            }
        }
    }
}

// Function to update particle velocities and positions
// This function implements the core PSO update rules:
// - Velocity is influenced by personal and global bests.
// - Position is updated based on velocity.
// - Boundary conditions are handled by reflecting particles when they exceed bounds.
void update_particles(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds,
                      ObjectiveFunction objective_function, double *global_best_position, double *global_best_value, double w, double c1, double c2) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            double r1 = random_double(0.0, 1.0); // Random coefficient for personal best influence
            double r2 = random_double(0.0, 1.0); // Random coefficient for global best influence

            // Update velocity using the PSO equation
            particles[i].velocity[j] = w * particles[i].velocity[j]
                                     + c1 * r1 * (particles[i].best_position[j] - particles[i].position[j])
                                     + c2 * r2 * (global_best_position[j] - particles[i].position[j]);

            // Clamp velocity to a max range to prevent wild swings
            double max_velocity = (bounds[j].upperBound - bounds[j].lowerBound) * 0.2; // 20% of the search space range
            if (particles[i].velocity[j] > max_velocity) particles[i].velocity[j] = max_velocity;
            if (particles[i].velocity[j] < -max_velocity) particles[i].velocity[j] = -max_velocity;

            // Update position
            particles[i].position[j] += particles[i].velocity[j];

            // Reflect position if it goes out of bounds
            if (particles[i].position[j] < bounds[j].lowerBound) {
                particles[i].position[j] = bounds[j].lowerBound + (bounds[j].lowerBound - particles[i].position[j]);
                particles[i].velocity[j] *= -1; // Reverse velocity direction
            }
            if (particles[i].position[j] > bounds[j].upperBound) {
                particles[i].position[j] = bounds[j].upperBound - (particles[i].position[j] - bounds[j].upperBound);
                particles[i].velocity[j] *= -1; // Reverse velocity direction
            }
        }

        // Evaluate fitness
        particles[i].value = objective_function(NUM_VARIABLES, particles[i].position);

        // Update personal best
        if (particles[i].value < particles[i].best_value) {
            particles[i].best_value = particles[i].value;
            for (int j = 0; j < NUM_VARIABLES; j++) {
                particles[i].best_position[j] = particles[i].position[j];
            }
        }

        // Update global best
        if (particles[i].value < *global_best_value) {
            *global_best_value = particles[i].value;
            for (int j = 0; j < NUM_VARIABLES; j++) {
                global_best_position[j] = particles[i].position[j];
            }
        }
    }
}

// Free memory allocated for particles
// This is a utility function to clean up dynamically allocated memory for particle attributes.
// Without this, the program would have memory leaks.
void free_particles(Particle *particles, int NUM_PARTICLES) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        free(particles[i].position);
        free(particles[i].velocity);
        free(particles[i].best_position);
    }
}

// PSO implementation
// This function orchestrates the entire Particle Swarm Optimization process:
// - Particles are initialized with random positions and velocities.
// - Iteratively, particles are updated to explore the search space.
// - The best solution found by the swarm is returned.
double pso(ObjectiveFunction objective_function, int NUM_VARIABLES, Bound *bounds, int NUM_PARTICLES, int MAX_ITERATIONS, double *best_position) {
    // Allocate memory for particles and global best variables
    Particle *particles = malloc(NUM_PARTICLES * sizeof(Particle));
    double global_best_value = DBL_MAX; // Start with the worst possible global best value
    double *global_best_position = malloc(NUM_VARIABLES * sizeof(double));

    // Initialize particles
    initialize_particles(particles, NUM_PARTICLES, NUM_VARIABLES, bounds, objective_function, global_best_position, &global_best_value);

    // Main PSO loop
    double w = 0.7, c1 = 1.5, c2 = 1.5; // PSO hyperparameters
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        update_particles(particles, NUM_PARTICLES, NUM_VARIABLES, bounds, objective_function,
                         global_best_position, &global_best_value, w, c1, c2);
    }

    // Copy global best position to output
    for (int i = 0; i < NUM_VARIABLES; i++) {
        best_position[i] = global_best_position[i];
    }

    // Free memory
    free_particles(particles, NUM_PARTICLES);
    free(particles);
    free(global_best_position);

    return global_best_value;
}
