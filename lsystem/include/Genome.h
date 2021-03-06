//
// Created by Karine Miras on 07/03/2017.
//

#ifndef LSYSTEM_PROTO_GENOME_H
#define LSYSTEM_PROTO_GENOME_H

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QGraphicsView>
#include <QtGui/QApplication>
#include <QtGui/QGraphicsRectItem>
#include <QtGui/QGraphicsLineItem>



#include "DecodedGeneticString.h"
#include "GeneticString.h"
#include "LSystem.h"



/**
 * Genome of an individual.
 */

class Genome{

public:

    Genome(std::string _id,
           std::string _id_parent1,
           std::string _id_parent2){
        id = _id;
        id_parent1 = _id_parent1;
        id_parent2 = _id_parent2;
        initializer();
    }


    unsigned int getTo();

    void  initializer();

    GeneticString  build_genetic_string(GeneticString  &gs,
                                         std::vector<std::string> genetic_string_items);

    GeneticString  getGeneticString();

    void setGeneticString(GeneticString  _gs);

    std::vector<std::string> getAxiom();
    void generate_final_string(int  replacement_iterations,
                               int export_genomes,
                               int generation,
                               std::string path,
                               int post_measuring);

    void decodeGeneticString(LSystem LS,
                             std::map<std::string, double> params,
                             std::string path);

    void constructor(int argc,
                     char* argv[],
                     std::map<std::string,
                             double> params,
                     std::string path,
                     int post_measuring);

    void draw_component( std::string parent_convertion,
                        int convertion_level,
                        std::string _directoryPath,
                        std::string reference,
                        std::string direction,
                        QGraphicsScene * scene,
                        std::vector<QGraphicsRectItem *>  items,
                        DecodedGeneticString::Vertex * c1,
                        DecodedGeneticString::Vertex * c2,
                        std::map<std::string, double> params);
    std::string getId();

    std::string getId_parent1();

    std::string getId_parent2();

    double getFit_parent1();

    double getFit_parent2();

    int getValid();

    void exportGenome(std::string path);

    void createEmbryo();

    void developGenomeIndirect(int argc,
                               char* argv[],
                               std::map<std::string, double> params,
                               LSystem LS,
                               int generation,
                               std::string path,
                               int post_measuring);

    void developGenomeDirect(int argc,
                             char* argv[],
                             std::map<std::string, double> params,
                             LSystem LS,
                             int generation,
                             std::string path);

    DecodedGeneticString getDgs();

    std::map< std::string, double > getMeasures();

    std::map< std::pair<int, int>, std::string >  getList_components();

    void updateMeasure(std::string key, double value);

    double getBalanceFitness();

    double getLocomotionFitness();

    double getLocomotionReal();

    double getNoveltyFitness();

    double getPenaltyFitness();

    double getNoveltyLocomotionFitness();

    double getFinalFitness();

    std::map< std::string, GeneticString  > getGrammar();

    void setGrammar(std::map< std::string, GeneticString > grammar);

    void removeMeasure(std::string key);

    void updateBalanceFitness(double fitness);

    void updateLocomotionFitness(double x, double y, double z, double time, std::map< std::string, double > params);

    void updateLocomotionReal(double value);

    void updateNoveltyFitness(double fitness);

    void updatePenaltyFitness(double fitness);

    void updateNoveltyLocomotionFitness(double fitness);

    void updateFinalFitness(double fitness);


    void addLetterGrammar(std::string letter,
                          GeneticString  lgs);

    void convertYamlBody(   std::string parent_convertion,
                            std::string dirpath,
                            int convertion_level,
                            std::string direction,
                            DecodedGeneticString::Vertex * c2,
                            std::string sensor);

    void convertYamlBrain(std::string _directoryPath);

    void build_grammar(LSystem LS, std::map<std::string, double> params);

    void build_genome_direct(LSystem LS, std::map<std::string, double> params);

    void setValid(int valid);

    QApplication * getApp();

    QGraphicsScene * getScene();



protected:

    std::string id; // id identifying the genome

    std::string id_parent1; // id of parent1 of genome

    std::string id_parent2; // id of parent2 of genome

    double locomotion_fitness = 0 ;

    double locomotion_real = 0 ;

    double balance_fitness = 0;

    double novelty_fitness = 0;

    double novelty_locomotion_fitness = 0;

    double penalty_fitness = 0;

    double final_fitness = 0;


    int valid = 1; // valid 1=yes, 0=no

    GeneticString  gs =  GeneticString();

    DecodedGeneticString dgs; // graph that logically represents the connections among the components forming the body

    QApplication * app = NULL;

    QGraphicsScene * scene = NULL; // scene holding the phenotype

    // letter(s) of the alphabet to build the initial developmental state of the genetic-string
    std::vector<std::string> axiom;

    // genetic-strings of production rules for each letter of the alphabet
    std::map< std::string, GeneticString  > grammar;

    // list of body metrics about the genome
    std::map< std::string, double > measures;

    // list of all components of the body, keys are coordinates, value is a letter
    std::map< std::pair<int, int>, std::string > list_components;

};

#endif //LSYSTEM_PROTO_GENOME_H