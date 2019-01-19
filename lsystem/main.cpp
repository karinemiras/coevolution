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

    std::string type_script = "post_measuring";

    // Experiments with random search (no simulation)
    //  ps: simulation experiments (with fitness) are started from revolve, include this as a library
    if(type_script == "random_search")
    {
        EvolutionIndirect evolve_generation = EvolutionIndirect("random_1", "../../");

        evolve_generation.setupEvolution(); // LATER, read numrums directly from configuration file!

        for (int e = 1; e <= evolve_generation.getParams()["num_runs"]; e++) {

            EvolutionIndirect evolve_generation2 = EvolutionIndirect("random_" + std::to_string(e), "../../");
            // LATER check if folder exists before and load params
            evolve_generation2.setupEvolution();
            int load_generation = 0;
            int ini = 1; // LATER read from file for load!!

            for (int i = ini; i <= evolve_generation2.getParams()["num_generations"]; i++) {
                evolve_generation2.runExperiment_part1(i, load_generation);
                evolve_generation2.runExperiment_part2(i);
                load_generation = 0;
            }
        }
    }

    // Post processing of measures for old experiments
    if(type_script == "post_measuring")
    {
        std::string path_exps = "paper4/plain_s1_old/plain_s1";

        EvolutionIndirect evolve_generation = EvolutionIndirect("", "../../");
        evolve_generation.readParams();

        for (int e = 4; e <= evolve_generation.getParams()["num_runs"]; e++)
        {
            EvolutionIndirect evolve_generation2 = EvolutionIndirect(path_exps + "_" + std::to_string(e), "../../");
            evolve_generation2.postMeasuring();

        }
    }


    return 0;
}


