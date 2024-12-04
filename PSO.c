#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h> // For printf
#include "utility.h"

// Helper function to generate random numbers in a range
double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Function to initialize particles
void initialize_particles(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds, ObjectiveFunction objective_function, double *global_best_position, double *global_best_value) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].position = malloc(NUM_VARIABLES * sizeof(double));
        particles[i].velocity = malloc(NUM_VARIABLES * sizeof(double));
        particles[i].best_position = malloc(NUM_VARIABLES * sizeof(double));
        particles[i].best_value = DBL_MAX;

        for (int j = 0; j < NUM_VARIABLES; j++) {
            particles[i].position[j] = random_double(bounds[j].lowerBound, bounds[j].upperBound);
            particles[i].velocity[j] = random_double(-1.0, 1.0);
            particles[i].best_position[j] = particles[i].position[j];
        }

        particles[i].value = objective_function(NUM_VARIABLES, particles[i].position);
        if (particles[i].value < *global_best_value) {
            *global_best_value = particles[i].value;
            for (int j = 0; j < NUM_VARIABLES; j++) {
                global_best_position[j] = particles[i].position[j];
            }
        }
    }
}

// Function to update particle velocities and positions
void update_particles(Particle *particles, int NUM_PARTICLES, int NUM_VARIABLES, Bound *bounds,
                      ObjectiveFunction objective_function, double *global_best_position, double *global_best_value, double w, double c1, double c2) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            double r1 = random_double(0.0, 1.0);
            double r2 = random_double(0.0, 1.0);

            // Update velocity
            particles[i].velocity[j] = w * particles[i].velocity[j]
                                     + c1 * r1 * (particles[i].best_position[j] - particles[i].position[j])
                                     + c2 * r2 * (global_best_position[j] - particles[i].position[j]);

            // Clamp velocity
            double max_velocity = (bounds[j].upperBound - bounds[j].lowerBound) * 0.2;
            if (particles[i].velocity[j] > max_velocity) particles[i].velocity[j] = max_velocity;
            if (particles[i].velocity[j] < -max_velocity) particles[i].velocity[j] = -max_velocity;

            // Update position
            particles[i].position[j] += particles[i].velocity[j];

            // Reflect position if out of bounds
            if (particles[i].position[j] < bounds[j].lowerBound) {
                particles[i].position[j] = bounds[j].lowerBound + (bounds[j].lowerBound - particles[i].position[j]);
                particles[i].velocity[j] *= -1;  // Reverse direction
            }
            if (particles[i].position[j] > bounds[j].upperBound) {
                particles[i].position[j] = bounds[j].upperBound - (particles[i].position[j] - bounds[j].upperBound);
                particles[i].velocity[j] *= -1;  // Reverse direction
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
void free_particles(Particle *particles, int NUM_PARTICLES) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        free(particles[i].position);
        free(particles[i].velocity);
        free(particles[i].best_position);
    }
}

// PSO implementation
// This version includes an improved threshold mechanism for early stopping
double pso(ObjectiveFunction objective_function, int NUM_VARIABLES, Bound *bounds, int NUM_PARTICLES, int MAX_ITERATIONS, double *best_position) {
    // Allocate memory for particles and global best variables
    Particle *particles = malloc(NUM_PARTICLES * sizeof(Particle));
    double global_best_value = DBL_MAX;
    double *global_best_position = malloc(NUM_VARIABLES * sizeof(double));

    // Initialize particles
    initialize_particles(particles, NUM_PARTICLES, NUM_VARIABLES, bounds, objective_function, global_best_position, &global_best_value);

    // Main PSO loop
    double w = 0.7, c1 = 1.5, c2 = 1.5; // PSO hyperparameters
    double threshold = 1e-13;            // Threshold for early stopping
    int patience = 50;                  // Iterations to wait for improvement
    int convergence_counter = 0;        // Counter for consecutive non-improvement iterations

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        update_particles(particles, NUM_PARTICLES, NUM_VARIABLES, bounds, objective_function,
                         global_best_position, &global_best_value, w, c1, c2);

        // Early stopping condition
        if (global_best_value < threshold) {
            convergence_counter++;
            if (convergence_counter >= patience) {
                printf("Stopping early after %d iterations, fitness: %lf\n", iter + 1, global_best_value);
                break;
            }
        } else {
            // Reset counter if improvement is observed
            convergence_counter = 0;
        }
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
