//
// Created by Karine Miras on 21/03/2017.
//

#ifndef LSYSTEM_PROTO_EVOLUTION_H
#define LSYSTEM_PROTO_EVOLUTION_H

#include <map>
#include <string>
#include <vector>

#include "Aux.h"
#include "Genome.h"
#include "Measures.h"
#include "Tests.h"


/**
 * Evolutionary algorithms.
 */
class  Evolution{


public:

    explicit Evolution(std::string experiment_name,
                       std::string path){

        this->experiment_name = experiment_name;
        this->path = path;
        initializer();

        this->measures_names["branching"] = "branching";
        this->measures_names["connectivity1"] = "connectivity1";
        this->measures_names["connectivity2"] = "connectivity2";
        this->measures_names["coverage"] = "coverage";
        this->measures_names["effective_joints"] = "effective_joints";
        this->measures_names["length_ratio"] = "length_ratio";
        this->measures_names["sensors"] = "sensors";
        this->measures_names["symmetry"] = "symmetry";
        this->measures_names["total_components"] = "total_components";
        this->measures_names["amplitude_average"] = "amplitude_average";
        this->measures_names["amplitude_deviation"] = "amplitude_deviation";
        this->measures_names["offset_average"] = "offset_average";
        this->measures_names["offset_deviation"] = "offset_deviation";
        this->measures_names["period_average"] = "period_average";
        this->measures_names["period_deviation"] = "period_deviation";
        this->measures_names["intra_params_dev_average"] = "intra_params_dev_average";
        this->measures_names["inter_params_dev_average"] = "inter_params_dev_average";
        this->measures_names["inputs_reach"] = "inputs_reach";
        this->measures_names["recurrence"] = "recurrence";
        this->measures_names["synaptic_reception"] = "synaptic_reception";
        this->measures_names["total_components_abs"] = "total_components_abs";
        this->measures_names["joints_abs"] = "joints_abs";
        this->measures_names["connectivity4"] = "connectivity4";
        this->measures_names["connectivity1_abs"] = "connectivity1_abs";
        this->measures_names["connectivity2_abs"] = "connectivity2_abs";
        this->measures_names["connectivity3"] = "connectivity3";

    }

    void saveHistory(int generation);

    void readParams();
    void testGeneticString(int argc,
                           char* argv[],
                           std::string test_genome);
    void measureIndividuals(int generation,
                            std::vector<Genome>  &individuals,
                            std::string dirpath,
                            int post_measuring);


    void evaluateLocomotion(int generation,
                     std::vector<Genome >  &individuals);
    void savesValidity(int generation);

    int  tournament();
    int  tournament_single();
    int  tournament_multi();
    void selection();
    std::vector<Genome>  getPopulation();
    std::map<std::string, double> getParams();
    double runExperiment_part1(int generation, int load_experiment);
    double runExperiment_part2(int generation);
    void exportGenerationMetrics(int generation,
                                 std::vector<int> metrics);
    void saveLocomotionInfo(std::string genome_id,
                             double x, double y, double z, double time);
    void saveBalance(std::string genome_id, double fitness);
    void exportPop(int generation);
    void calculateNovelty();
    void postMeasuring();
    void calculateNoveltyLocomotion();
    void calculateFinalFitness();
    void calculatePenaltyFitness();
    void calculateRankFitness();
    void saveParameters();
    void surv_selection_tournament();
    void surv_selection_comma();
    void logsTime(std::string moment);
    void setupEvolution();
    void writesEvolutionState(int generation, int next_id);
    std::vector<std::string> readsEvolutionState();
    void loadsParams();
    void loadIndividuals(int generation, std::string type);
    std::vector<int> calculateNicheCoverage();
    void createHeader();
    void updateParameter(std::string key, double value);
    void developIndividuals(int argc, char* argv[],
                            LSystem LS,
                            int generation,
                            std::vector<Genome>  &individuals,
                            std::string dir,
                            int post_measuring);
    int loadExperiment();
    void initExperiment(int argc, char* argv[],
                       LSystem LS);
    void summaryNicheCoverage();
    int getGeneration_genome(std::string idgenome);
    double compareIndividual(Measures m,
                             std::string idgenome);
    double compareParents(std::string idparent1,
                          std::string idparent2);
    void  addToArchive();
    void cleanMemory(std::vector< int > index_selected);
    void cleanVertex(DecodedGeneticString::Vertex * v);

    virtual void initPopulation(LSystem LS){};
    virtual void crossover(LSystem LS){};
    virtual void mutation(LSystem LS){};
    void initializer();



protected:

    std::string experiment_name = ""; // name for the experiment

    std::string path = ""; // path of the lsystem

    std::map<std::string, double> params; // contains the list of parameters loaded from parameter file

    // containsgeneral auxiliar methods for the experiments
    Aux aux = Aux(this->experiment_name,this->getParams(),this->path);
    // contains methods with tests for the system
    Tests tests = Tests(this->experiment_name,
                        this->getParams(),
                        this->path);

    int next_id = 0; // id that will be given for the next genome to be created


    std::map< std::string, std::string > measures_names;

    // points in a grid representing the morphological space
    std::map<std::string, std::vector<double>> morphological_grid_generation;

    std::map<std::string, std::vector<std::string>> morphological_grid_accumulated;

    // contains the genomes of all the individuals of the current population
    std::vector<Genome>  population;

    // contains the genomes of all the individuals the new offspring
    std::vector<Genome>  offspring;

    // contains the genomes of all individuals in the archive
    std::vector<Genome>  archive;


};


#endif //LSYSTEM_PROTO_EVOLUTION_H
