#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "utility.h"
#include "OF_lib.h"

int main(int argc, char **argv) {

    // Validate arguments
    if (argc < 7) {
        fprintf(stderr, "Usage: %s <ObjectiveFunctionName> <NUM_VARIABLES> <LowerBound> <UpperBound> <NUM_PARTICLES> <MAX_ITERATIONS>\n", argv[0]);
        return 1;
    }

    // Parse command-line arguments
    char *objective_function_name = argv[1];
    int NUM_VARIABLES = atoi(argv[2]);
    double lower_bound = atof(argv[3]);
    double upper_bound = atof(argv[4]);
    int NUM_PARTICLES = atoi(argv[5]);
    int MAX_ITERATIONS = atoi(argv[6]);

    // Assign the appropriate function pointer based on user input
    ObjectiveFunction objective_function = NULL;

    if (strcmp(objective_function_name, "griewank") == 0) {
        objective_function = griewank;
    } else if (strcmp(objective_function_name, "levy") == 0) {
        objective_function = levy;
    } else if (strcmp(objective_function_name, "rastrigin") == 0) {
        objective_function = rastrigin;
    } else if (strcmp(objective_function_name, "rosenbrock") == 0) {
        objective_function = rosenbrock;
    } else if (strcmp(objective_function_name, "schwefel") == 0) {
        objective_function = schwefel;
    } else if (strcmp(objective_function_name, "dixon_price") == 0) {
        objective_function = dixon_price;
    } else if (strcmp(objective_function_name, "michalewicz") == 0) {
        objective_function = michalewicz;
    } else if (strcmp(objective_function_name, "styblinski_tang") == 0) {
        objective_function = styblinski_tang;
    } else {
        fprintf(stderr, "Invalid objective function: %s\n", objective_function_name);
        return 1;
    }

    // Print inputs for verification
    printf("Objective Function: %s\n", objective_function_name);
    printf("Number of Variables: %d\n", NUM_VARIABLES);
    printf("Lower Bound for all variables: %lf\n", lower_bound);
    printf("Upper Bound for all variables: %lf\n", upper_bound);

    // Allocate bounds
    Bound *bounds = (Bound *)malloc(NUM_VARIABLES * sizeof(Bound));
    for (int i = 0; i < NUM_VARIABLES; i++) {
        bounds[i].lowerBound = lower_bound;
        bounds[i].upperBound = upper_bound;
    }

    // Allocate memory for best position
    double *best_position = (double *)malloc(NUM_VARIABLES * sizeof(double));

    // Measure CPU time for PSO
    clock_t start_time = clock();
    double best_fitness = pso(objective_function, NUM_VARIABLES, bounds, NUM_PARTICLES, MAX_ITERATIONS, best_position);
    clock_t end_time = clock();
    double elapsed_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    // Print results
    printf("Optimal fitness: %lf\n", best_fitness);
    printf("Optimal position: ");
    for (int i = 0; i < NUM_VARIABLES; i++) {
        printf("%lf ", best_position[i]);
    }
    printf("\nCPU Time: %lf seconds\n", elapsed_time);

    // Free memory
    free(bounds);
    free(best_position);

    return 0;
}
