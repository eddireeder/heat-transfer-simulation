#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>


// Function to take a temperature map as input and generate a gnuplot script to display the data in a heat map
int output_heatmap(int square_array_length, float temperature_map[square_array_length][square_array_length]) {

  FILE *fp;
  fp = fopen("heatmap.gnu", "w");
  fprintf(fp, "set ylabel 'y (m)'\n");
  fprintf(fp, "set xlabel 'x (m)'\n");
  // A messy way to format the scales, but this isnt about how to correctly use gnuplot
  fprintf(fp, "set xtics ('0' 0, '0.02' 100, '0.04' 200, '0.06' 300, '0.08' 400, '0.1' 499)\n");
  fprintf(fp, "set ytics ('0' 0, '0.02' 100, '0.04' 200, '0.06' 300, '0.08' 400, '0.1' 499)\n");
  fprintf(fp, "$map2 << EOD\n");
  // Loop through the temperature map, and write the coordinates and heat value
  for (int i = 0; i < square_array_length; i++) {
    for (int j = 0; j < square_array_length; j++) {
      fprintf(fp, "%d %d %f\n", i, j, temperature_map[i][j]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "EOD\n");
  fprintf(fp, "plot '$map2' using 2:1:3 with image\n");
  fprintf(fp, "pause -1 'Hit any key to continue'");
  fclose(fp);

  return 1;
}


// Simulate the heat transfer equation and return the time taken to complete if using OpenMP
float simulate_hte(const int steps) {

  // Start timer
  #ifdef _OPENMP
    float tstart = omp_get_wtime();
  #endif
  
  // Initial conditions
  const float alpha_squared = 0.01;
  const int num_intervals = 500;
  const float domain_square_length = 0.1;
  const float hot_region_square_length = 0.04;
  const float background_temp = 20.0;
  const float hot_region_temp = 50.0;
  const float time_step = 0.0000009;
  const int num_threads = 4;

  // Initialise temperature map
  float temperature_map[num_intervals][num_intervals];
  const int hot_region_array_length = num_intervals*(hot_region_square_length/domain_square_length);
  const int hot_region_lower_bound = (num_intervals/2) - (hot_region_array_length/2);
  const int hot_region_upper_bound = (num_intervals/2) + (hot_region_array_length/2);
  #pragma omp parallel for default(none) shared(temperature_map)
    for (int i = 0; i < num_intervals; i++) {
      for (int j = 0; j < num_intervals; j++) {
        if (i >= hot_region_lower_bound && i <= hot_region_upper_bound & j >= hot_region_lower_bound && j <= hot_region_upper_bound) {
          temperature_map[i][j] = hot_region_temp;
        } else {
          temperature_map[i][j] = background_temp;
        }
      }
    }

  // Calculate distance between intervals
  const float delta_x = domain_square_length/num_intervals;
  const float delta_y = domain_square_length/num_intervals;

  // Compute time steps
  for (int i = 0; i < steps; i++) {
    printf("Computing step %d\n", i);
    // Create new 2D array to store updated values without affecting later calculations (stay in the same time step)
    float new_temperature_map[num_intervals][num_intervals];
    #pragma omp parallel for default(none) shared(temperature_map, new_temperature_map)
      for (int i = 0; i < num_intervals; i++) {
        for (int j = 0; j < num_intervals; j++) {
          // Assign surrounding temps, checking whether on edge of map (set temp so that no change comes from outside the map)
          float left_temp, right_temp, top_temp, bottom_temp;
          if (i == 0) {
            left_temp = temperature_map[i][j];
          } else {
            left_temp = temperature_map[i - 1][j];
          }
          if (i == (num_intervals - 1)) {
            right_temp = temperature_map[i][j];
          } else {
            right_temp = temperature_map[i + 1][j];
          }
          if (j == 0) {
            top_temp = temperature_map[i][j];
          } else {
            top_temp = temperature_map[i][j - 1];
          }
          if (j == (num_intervals - 1)) {
            bottom_temp = temperature_map[i][j];
          } else {
            bottom_temp = temperature_map[i][j + 1];
          }
          // Calculate the rate of change of temperature
          float horizontal_rate_of_change_of_temp = alpha_squared*((left_temp + right_temp - 2*temperature_map[i][j])/(delta_x*delta_x));
          float vertical_rate_of_change_of_temp = alpha_squared*((top_temp + bottom_temp - 2*temperature_map[i][j])/(delta_y*delta_y));
          float rate_of_change_of_temp = horizontal_rate_of_change_of_temp + vertical_rate_of_change_of_temp;
          // Calculate new temperature
          new_temperature_map[i][j] = temperature_map[i][j] + rate_of_change_of_temp*time_step;
        }
      }
    // Copy new temperature map as current
    memcpy(temperature_map, new_temperature_map, sizeof temperature_map);
  }

  printf("Done\n");

  // End timer
  #ifdef _OPENMP
    float time_to_complete = omp_get_wtime() - tstart;
    printf("Time taken to complete: %f\n", time_to_complete);    
  #endif

  // Write a gnuplot file
  output_heatmap(num_intervals, temperature_map);

  #ifdef _OPENMP
    return time_to_complete;
  #else
    return 0.0;
  #endif
}


// Main function
int main(int argc, char *argv[]) {

    // Initialise
    int steps = 10000;
    int num_threads = 1;
    bool investigate = false;

    // Loop through arguments and check for flags
    for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "--steps") == 0) {
        steps = atoi(argv[i + 1]);
      } else if (strcmp(argv[i], "--num-threads") == 0) {
        num_threads = atoi(argv[i + 1]);
      } else if (strcmp(argv[i], "--investigate") == 0) {
        investigate = true;
      }
    }

    if (investigate) {
      // Run simulation multiple times and produce graph
      #ifdef _OPENMP
        printf("Starting investigation\n");
        printf("There are %d processors available\n", omp_get_num_procs());
        const int max_threads = omp_get_num_procs();
        float results[max_threads];
        for (int i = 0; i < max_threads; i++) {
          printf("Setting number of threads to %d\n", (i + 1));
          // Set number of threads to use
          omp_set_num_threads(i + 1);
          results[i] = simulate_hte(steps);
        }
        for (int i = 0; i < max_threads; i++) {
          printf("%f\n", results[i]);
        }
      #else
        printf("Not using OpenMP, can't investigate parallel performance\n");
      #endif

      return 1;
    }

    // Ensure program compiles and runs with or without OpenMP
    #ifdef _OPENMP
      printf("Using OpenMP\n");
      printf("There are %d processors available\n", omp_get_num_procs());
      printf("Setting number of threads to %d\n", num_threads);
      // Set number of threads to use
      omp_set_num_threads(num_threads);
    #else
      printf("Not using OpenMP\n");
    #endif

    // Run
    simulate_hte(steps);
    return 1;
}