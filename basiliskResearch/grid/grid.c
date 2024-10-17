#include "grid/multigrid.h"
#include "run.h"

int main() {
  size (1.0);           // Sets the domain size
  init_grid (64);       // Initializes the grid
  run();                // Runs the simulation

  return 0;
}
