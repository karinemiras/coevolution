//
// Created by Karine Miras on 21/03/2017.
//

#ifndef LSYSTEM_PROTO_EVOLUTIONINDIRECT_H
#define LSYSTEM_PROTO_EVOLUTIONINDIRECT_H

#include <map>
#include <string>
#include <vector>

#include "Aux.h"
#include "Evolution.h"
#include "Genome.h"
#include "Tests.h"

/**
 * Evolutionary algorithm.
 */

class EvolutionIndirect: public Evolution{


public:


    EvolutionIndirect(std::string experiment_name,
                      std::string path) :
            Evolution(experiment_name,
                      path){ }


    void initPopulation(LSystem LS);
    void crossover(LSystem LS);
    void mutation(LSystem LS);



};


#endif //LSYSTEM_PROTO_EVOLUTIONINDIRECT_H
