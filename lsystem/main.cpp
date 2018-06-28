/*
 * Genome indirect enconding: L-System
 * Author: Karine Miras
 * Created: 02/03/17
 */

#include<iostream>
#include <random>
#include "EvolutionIndirect.h"

using namespace std;


int main(int argc, char* argv[])
{
    EvolutionIndirect evolve_generation = EvolutionIndirect("test","../../");

    /* test setup */
      evolve_generation.setupEvolution();
      int load_generation = 0;
      int ini = 1;
    /* test setup */

    for(int i=ini; i <= evolve_generation.getParams()["num_generations"]; i++)
    {
        evolve_generation.runExperiment_part1(i, load_generation);
        evolve_generation.runExperiment_part2(i);
        load_generation = 0;
    }

    return 0;
}


