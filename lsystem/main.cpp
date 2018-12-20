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


    EvolutionIndirect evolve_generation = EvolutionIndirect("random_1","../../");

    evolve_generation.setupEvolution();

    for(int e=1; e <= evolve_generation.getParams()["num_runs"]; e++)
    {
    
          evolve_generation = EvolutionIndirect("random_"+std::to_string(e),"../../");

          evolve_generation.setupEvolution();
          int load_generation = 0;
          int ini = 1; // read from file later for load!!

        for(int i=ini; i <= evolve_generation.getParams()["num_generations"]; i++)
        {
            evolve_generation.runExperiment_part1(i, load_generation);
            evolve_generation.runExperiment_part2(i);
            load_generation = 0;
        }
    }


    return 0;
}


