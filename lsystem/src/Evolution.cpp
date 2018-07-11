//
// Created by Karine Miras on 21/03/2017.
//

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>


#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>

using namespace mlpack;
using namespace mlpack::neighbor; // NeighborSearch and NearestNeighborSort
using namespace mlpack::metric; // EuclideanDistance

#include "Aux.h"
#include "Evolution.h"
#include "Genome.h"
#include "LSystem.h"
#include "Measures.h"


/**
 * Reads parameters from file.
 **/
void Evolution::readParams()
{
  std::string line;
  std::ifstream myfile(this->path+"lsystem/configuration.txt");

  /*   pop_size - size of the population of genomes
  *    offspring_prop - proportion of the population size to dcalculate size of offspring
  *    num_initial_comp - number of initial (random) components in the production rules of the grammar
  *    show_phenotypes - flag to show the phenotype graphic
  *    export_phenotypes - if exports the phenotypes to images (1) or not (0)
  *    export_genomes - if exports the genomes to files (1) or not (0)
  *    replacement_iterations - number of replacement iterations for the l-system
  *    size_component - size of each component in pixels
  *    spacing - spacing between components in pixels
  *    mutation_prob - probability of adding/removing/swaping items (letters/commands) to the genetic-string in the mutation
  *    max_comps - maximum number of components allowed per phenotype
  *    prob_add_archive - probability of adding any genome to the archive
  *    grid_bins - number of bins to break morphological space
  *    logs_to_screen - if exports the logs to the screen (1) or not (0)
  *    logs_to_file - if exports logs to a file (1) or not (0)
  *    max_attempt_multi  - maximum number of tournaments for multiobjetive
  *    objective - 0 for single and 1 for multi
  */

  if (myfile.is_open())
  {
    while (getline(
        myfile,
        line))
    {
      std::vector< std::string > tokens;

      // parameters label and value separated by space
      boost::split(
          tokens,
          line,
          boost::is_any_of(" "));

      // first item is the label, second is the value
      this->params[tokens[0]] = std::stod(tokens[1]);
    }
    myfile.close();
  }

}

/*
 * Change value of parameter for an experiment.
 * */
void Evolution::updateParameter(
    std::string key,
    double value)
{

  this->params[key] = value;
}


/**
 * Loads parameters from saved state.
 **/
void Evolution::loadsParams()
{

  std::string line;
  std::ifstream myfile(
      this->path+"experiments/"
      + this->experiment_name +
      "/configuration.txt");

  if (myfile.is_open())
  {
    while (getline(
        myfile,
        line))
    {

      std::vector< std::string > tokens;
      // parameters label and value separated by space
      boost::split(
          tokens,
          line,
          boost::is_any_of(" "));

      // first item is the label, second is the value
      this->params[tokens[0]] = std::stod(tokens[1]);
    }
    myfile.close();
  }
  else
  {
    this->aux.logs("Unable to open parameters state.");
  }

}



/* Finds out in which generation the genome was generated.
 * @param idgenome - id of the genome for which to verify generation.
 * */
int Evolution::getGeneration_genome(std::string idgenome)
{

  int generation_genome = 0;
  int offspring_size =
      this->params["pop_size"] * this->params["offspring_prop"];

  // generation of the genome can be found by its id, considering the size of the population and the offspring
  if (this->params["offspring_prop"] == 1)
  {
    generation_genome = (int) trunc(
        std::stoi(idgenome)
        / this->params["pop_size"]) + 1;
  }
  else
  {
    generation_genome = (int) trunc((std::stoi(idgenome) - offspring_size)
                                    / offspring_size) + 1;
  }

  if (generation_genome == 0)
  { generation_genome = 1; }

  return generation_genome;

}


/**
*  Copies the phenotypes of the selection population to a separate folder.
*  @param generation - number of generation in evolution
**/
void Evolution::exportPop(int generation)
{
  std::ofstream measures_file_general;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/measures2.txt";
  measures_file_general.open(
      path,
      std::ofstream::app);

  for (int i = 0;
       i < this->population.size();
       i++)
  { // for each genome in the population


    // finds number of generation to which the genome belongs to
    int generation_genome = this->getGeneration_genome(this->population[i].getId());

    std::string filename = "/body_" + this->population[i].getId()
                           + "_p1_" +
                           this->population[i].getId_parent1() +
                           "_p2_" +
                           this->population[i].getId_parent2() +
                           ".png";

    std::string filename2 = "/brain_" + this->population[i].getId()+ ".png";

    std::string pathfrom = this->path+"experiments/"
                           + this->experiment_name + "/offspringpop" +
                           std::to_string(generation_genome);

    std::string pathto = this->path+"experiments/"
                         + this->experiment_name + "/selectedpop" +
                         std::to_string(generation);

    // copies phenotype file from offspring folder to selected population folder
    system(("exec cp " + pathfrom + filename + " " + pathto +
            filename).c_str());
    system(("exec cp " + pathfrom + filename2 + " " + pathto +
            filename2).c_str());


    // copies values of metrics to file of selected population
    std::string line;
    std::ifstream measures(
        this->path+"experiments/" + this->experiment_name +
        "/offspringpop" +
        std::to_string(generation_genome)
        + "/measures" +
        this->population[i].getId() + ".txt");

    while (getline(
        measures,
        line))
    {

      std::vector< std::string > tokens;
      boost::split(
          tokens,
          line,
          boost::is_any_of(":"));

      measures_file_general << std::to_string(generation) << " "
                            << this->population[i].getId() << " "
                            << tokens[0] << " " << tokens[1] << std::endl;
    }
  }

  measures_file_general.close();
}


/*
 * Compare phenotype of the individual with its parent's.
 * */
double Evolution::compareIndividual(
    Measures m,
    std::string idgenome)
{


  int generation_genome = this->getGeneration_genome(idgenome);

  std::string line;
  std::ifstream measures(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" + std::to_string(generation_genome) +
      "/measures" + idgenome + ".txt");

  double dif = 0;
  while (getline(
      measures,
      line))
  {

    std::vector< std::string > tokens;
    boost::split(
        tokens,
        line,
        boost::is_any_of(":"));

    dif += std::pow(
        m.getGenome()->getMeasures()[tokens[0]] - std::stod(tokens[1]),
        2);
  }
  dif = roundf(std::sqrt(dif) * 100) / 100;

  return dif;

}


/*
 * Compare phenotype of the parents.
 * */
double Evolution::compareParents(
    std::string idparent1,
    std::string idparent2)
{


  int generation_genome_parent1 = this->getGeneration_genome(idparent1);
  std::string line;
  std::ifstream measures(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" +
      std::to_string(generation_genome_parent1) +
      "/measures" + idparent1 + ".txt");

  int generation_genome_parent2 = this->getGeneration_genome(idparent2);
  std::string line2;
  std::ifstream measures2(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" +
      std::to_string(generation_genome_parent2) +
      "/measures" + idparent2 + ".txt");

  double dif = 0;
  while (getline(
      measures,
      line))
  {

    getline(
        measures2,
        line2);

    std::vector< std::string > tokens, tokens2;
    boost::split(
        tokens,
        line,
        boost::is_any_of(":"));
    boost::split(
        tokens2,
        line2,
        boost::is_any_of(":"));

    dif += std::pow(
        std::stod(tokens[1]) - std::stod(tokens2[1]),
        2);
  }
  dif = roundf(std::sqrt(dif) * 100) / 100;

  return dif;

}

/**
 * Measures all the individuals of the population for several metrics.
 *  @param argc - command line parameter
 *  @param argv[] - command line parameter
 *  @param individuals - array with genomes
 *  @param dirpath - name of the output directory
 **/
void Evolution::measureIndividuals(
    int generation,
    std::vector< Genome > &individuals,
    std::string dirpath)
{

  std::ofstream differences_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/distances.txt";
  differences_file.open(
      path,
      std::ofstream::app);


  // for each genome of the population
  for (int i = 0;
       i < individuals.size();
       i++)
  {
    std::cout<<"idgen"<<individuals[i].getId()<<std::endl;

    Measures m = Measures(
        this->experiment_name,
        this->params,
        this->path
    );
    m.setGenome(individuals[i]);
    // measures phenotype
    m.measurePhenotype(
        this->params,
        dirpath,
        generation);
    std::cout<<"mediu"<<std::endl;
    // compares measures between individuals
    if (individuals[i].getId_parent1() != "N")
    {

      double dif = this->compareIndividual(
          m,
          individuals[i].getId_parent1());
      differences_file << individuals[i].getId() << " " << dif;

      dif = this->compareIndividual(
          m,
          individuals[i].getId_parent2());
      differences_file << " " << dif;

      dif = this->compareParents(
          individuals[i].getId_parent1(),
          individuals[i].getId_parent2());
      differences_file << " " << dif << std::endl;

    }
    std::cout<<"comparou"<<std::endl;
  }

  differences_file.close();
}



/*
 *
 * */




/**
 * Creates files of results containing headers.
 */

void Evolution::createHeader()
{

  std::ofstream file;

  std::string path =
      this->path+"experiments/" + this->experiment_name + "/history.txt";
  file.open(path);
  file
      << "generation idgenome noveltyfit locofit finalfit balfit "
          "novlfit penfit idparent1 idparent2 "
      << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/evolution.txt";
  file.open(path);
  file
      << "generation idbest_nov maxfit_nov meanfit_nov idbest_loco "
          "maxfit_loco meanfit_loco idbest_fin maxfit_fin meanfit_fin "
          "idbest_bal maxfit_bal meanfit_bal "
          "idbest_novl maxfit_novl meanfit_novl "
          "idbest_pen maxfit_pen meanfit_pen "
          "nichecoverage_generation nichecoverage_accumulated";
  file << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/measures.txt";
  file.open(path);
  file << "generation idgenome";
  for (const auto &it : this->measures_names)
  {
    file << " " << it.first;
  }
  file << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/measures2.txt";
  file.open(path);
  file << "generation genome measures value" << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/distances.txt";
  file.open(path);
  file << "idgenome distances_parent1 distances_parent2 distances_parents"
       << std::endl;
  file.close();


}

void Evolution::saveHistory(int generation)
{

  for (int i = 0;
       i < this->population.size();
       i++)
  {

    std::ofstream history_file;
    std::string path =
        this->path + "experiments/" + this->experiment_name + "/history.txt";
    history_file.open(
        path,
        std::ofstream::app);

    history_file << std::to_string(generation) << " "     // generation
                 << this->population[i].getId() << " "   // idgenome
                 << this->population[i].getNoveltyFitness() << " "
                 << this->population[i].getLocomotionFitness() << " "
                 << this->population[i].getFinalFitness() << " "
                 << this->population[i].getBalanceFitness() << " "
                 << this->population[i].getNoveltyLocomotionFitness() << " "
                 << this->population[i].getPenaltyFitness() << " "
                 << this->population[i].getId_parent1() << " "  // id of parent1
                 << this->population[i].getId_parent2() << " " // id of parent2
                 << std::endl;

    history_file.close();
  }
}



void Evolution::savesValidity(int generation)
{
  // saves list of robots and its validits to file to be read by simulator
  std::ofstream file;
  std::string path =
      this->path + "experiments/" + this->experiment_name +
      "/offspringpop" + std::to_string(generation) + "/validity_list.txt";
  file.open(path);
  for (int i = 0;
       i <  this->offspring.size();
       i++)
  {
    file << this->offspring[i].getId() << " " <<  this->offspring[i].getValid()
         << std::endl;
  }
  file.close();
}


/**
 *
 * @return - the index of the winner genome
 */

int Evolution::tournament_single()
{

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution< int > dist_1(
      0,
      (int) this->population.size() -
      1); // size of current pop (parents+offspring)

  int genome1 = dist_1(generator); // random genome 1
  int genome2 = dist_1(generator); // random genome 2

  // return the genome with higher fitness
  if (this->population[genome1].getFinalFitness() >
      this->population[genome2].getFinalFitness())
  {

    return genome1;
  }
  else
  {
    return genome2;
  }

}

int Evolution::tournament()
{
  if (this->getParams()["objective"] == 0){
    return this->tournament_single();
  }
  else{
    return this->tournament_multi();
  }
}


int Evolution::tournament_multi()
{

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution< int > dist_1(0,
      (int) this->population.size() - 1); // size of current pop (parents+offspring)

  int genome1 = -1; // random genome 1
  int genome2 = -1; // random genome 2
  int dominant = -1;
  int attempts = 0;


  while( dominant == -1 and attempts <= this->getParams()["max_attempt_multi"] )
  {

    genome1 = dist_1(generator);
    genome2 = dist_1(generator);

    attempts = attempts + 1;

    // genome1 dominates
    if (this->population[genome1].getBalanceFitness()
          >= this->population[genome2].getBalanceFitness()
        and this->population[genome1].getLocomotionFitness()
          >= this->population[genome2].getLocomotionFitness()
        and this->population[genome1].getNoveltyFitness()
          >= this->population[genome2].getNoveltyFitness()
        and
          (
              this->population[genome1].getBalanceFitness()
                 > this->population[genome2].getBalanceFitness()
              or this->population[genome1].getLocomotionFitness()
                 > this->population[genome2].getLocomotionFitness()
              or this->population[genome1].getNoveltyFitness()
                 > this->population[genome2].getNoveltyFitness()
          )
        ) dominant = genome1;

    // genome2 dominates
    if (this->population[genome2].getBalanceFitness()
            >= this->population[genome1].getBalanceFitness()
        and this->population[genome2].getLocomotionFitness()
            >= this->population[genome1].getLocomotionFitness()
        and this->population[genome2].getNoveltyFitness()
            >= this->population[genome1].getNoveltyFitness()
        and
          (
              this->population[genome2].getBalanceFitness()
                 > this->population[genome1].getBalanceFitness()
              or this->population[genome2].getLocomotionFitness()
                 > this->population[genome1].getLocomotionFitness()
              or this->population[genome2].getNoveltyFitness()
                 > this->population[genome1].getNoveltyFitness()
          )
        ) dominant = genome2;
  }

  // if no dominent was found, chooses the first (random)
  if(dominant == -1) dominant = genome1;

  return dominant;

}

//
//int Evolution::tournament_multi()
//{
//
//  std::random_device rd;
//  std::default_random_engine generator(rd());
//  std::uniform_int_distribution< int > dist_1(0,
//                                              (int) this->population.size() - 1); // size of current pop (parents+offspring)
//
//  int genome1 = -1; // random genome 1
//  int genome2 = -1; // random genome 2
//  int dominant = -1;
//  int attempts = 0;
//
//  std::cout<<"-----------------"<<std::endl;
//  while( dominant == -1 and attempts <= this->getParams()["max_attempt_multi"] )
//  {
//
//    genome1 = dist_1(generator);
//    genome2 = dist_1(generator);
//
//    attempts = attempts + 1;
//
//    std::cout<<attempts<<std::endl;
//    std::cout<<this->population[genome1].getPenaltyFitness()<<std::endl;
//    std::cout<<this->population[genome2].getPenaltyFitness()<<std::endl;
//    std::cout<<this->population[genome1].getLocomotionFitness()<<std::endl;
//    std::cout<<this->population[genome2].getLocomotionFitness()<<std::endl;
//
//    // genome1 dominates
//    if (this->population[genome1].getPenaltyFitness()
//        > this->population[genome2].getPenaltyFitness()
//        and this->population[genome1].getLocomotionFitness()
//            > this->population[genome2].getLocomotionFitness()
//            and this->population[genome2].getNoveltyFitness()
//                > this->population[genome1].getNoveltyFitness()
//        ) dominant = genome1;
//
//    // genome2 dominates
//    if (this->population[genome2].getPenaltyFitness()
//        > this->population[genome1].getPenaltyFitness()
//        and this->population[genome2].getLocomotionFitness()
//            > this->population[genome1].getLocomotionFitness()
//            and this->population[genome2].getNoveltyFitness()
//                >= this->population[genome1].getNoveltyFitness()
//
//        ) dominant = genome2;
//  }
//
//  // if no dominent was foujd, chooses the first (random)
//  if(dominant == -1) dominant = genome1;
//
//  return dominant;
//
//}
//
//

/**
*  Selection of genomes in a population.
**/

void Evolution::selection()
{
  std::vector< Genome > selected = std::vector< Genome >();
  std::vector< int > index_selected = std::vector< int >();

  // selects genomes, maintaining population size
  for (int i = 0;
       i < this->params["pop_size"];
       i++)
  {

    int genome = -1;
    // selects one genome by tournament
    genome = this->tournament();


    // makes sure that the same genome wont be selected more than once
    while (std::find(
        index_selected.begin(),
        index_selected.end(),
        genome) != index_selected.end())
    {
      genome = this->tournament();
    }
    selected.push_back(this->population[genome]);
    index_selected.push_back(genome);
  }

  this->cleanMemory(index_selected);

  this->population = selected; // substitutes current population for the selected subset

  // # TEST: Tests if population size remains correct.
  this->tests.testPopsize(
      this->population,
      (int) this->params["pop_size"]);
}

/*
 * Deallocate memory used by the non-selected individuals.
 * */

void Evolution::cleanMemory(std::vector< int > index_selected)
{
  // cleaning memory
  std::cout<<"START cleaning memory for non-selected individuals"<<std::endl;
  for (int i = 0;
       i < this->population.size();
       i++)
  {
    // for non-selected individuals
    if(std::find(
        index_selected.begin(),
        index_selected.end(),
        i) == index_selected.end())
    {
      auto item = this->population[i].getGeneticString().getStart();
      while (item not_eq NULL)
      {
        auto item2 = item->next;
        delete item;
        item = item2;
      }

//      std::cout<<"grammar genetic-strings"<<std::endl;
//      for( auto &g: this->population[i].getGrammar())
//      {
//        item = g.second.getStart();
//        while (item not_eq NULL)
//        {
//          auto item2 = item->next;
//          delete item;
//          item = item2;
//        }
//      }

      this->cleanVertex(this->population[i].getDgs().getRoot());

      if(this->population[i].getScene() != NULL)
      {
        QList< QGraphicsItem * > all = this->population[i].getScene()->items();
        for (int i = 0; i < all.size(); i++)
        {
          QGraphicsItem *gi = all[i];
          if (gi->parentItem() == NULL) delete gi;
        }
        delete this->population[i].getScene();
      }
    }
  }
  std::cout<<"FINISH cleaning memory for non-selected individuals"<<std::endl;
}

void Evolution::cleanVertex(DecodedGeneticString::Vertex * v){

  if(v != NULL)
  {
    this->cleanVertex(v->left);
    this->cleanVertex(v->front);
    this->cleanVertex(v->right);
    if(v->item == "C")
      this->cleanVertex(v->back);
    delete v;
  }
}

/**
 *  Saves state of the generations to file.
 */

void Evolution::exportGenerationMetrics(
    int generation,std::vector<int> metrics)
{
  std::ofstream evolution_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/evolution.txt";
  evolution_file.open(
      path,
      std::ofstream::app);

  evolution_file << generation;

  // fetches all types of fitness
  for(int f=0; f<6; f++)
  {
    double maximum_fitness = 0;
    std::string best_genome = "0";
    double average_fitness = 0;
    double fitness = 0;

    for (int i = 0;
         i < this->getPopulation().size();
         i++)
    {
      if(f==0)
       fitness = this->getPopulation()[i].getNoveltyFitness();
      if(f==1)
        fitness = this->getPopulation()[i].getLocomotionFitness();
      if(f==2)
        fitness = this->getPopulation()[i].getFinalFitness();
      if(f==3)
        fitness = this->getPopulation()[i].getBalanceFitness();
      if(f==4)
        fitness = this->getPopulation()[i].getNoveltyLocomotionFitness();
      if(f==5)
        fitness = this->getPopulation()[i].getPenaltyFitness();

      // finds the maximum/best fitness of the population
      if (fitness > maximum_fitness)
      {
        best_genome = this->getPopulation()[i].getId();
        maximum_fitness = fitness;
      }
      average_fitness += fitness;
    }

    // calculates the average
    average_fitness /= this->getPopulation().size();

    evolution_file << " " << best_genome
                   << " " << maximum_fitness
                   << " " << average_fitness;
  }

  for (const auto &m : metrics)
  {
    evolution_file << " " << m;
  }

  evolution_file << std::endl;
  evolution_file.close();
}


void Evolution::setupEvolution()
{

  this->readParams();

  // cleans old files and creates folders for the experiment
  aux.removeFolder(this->path+"experiments/"+this->experiment_name);
  aux.createFolder(this->path+"experiments/"+this->experiment_name);

  // logs parameters configuration
  this->saveParameters();

}


/**
*  Develops genomes of the population: 1- grows genetic-string of the genome according to grammar, 2- decodes it, 3- constructs the phenotype
*  @param argc - default argument
*  @param argv[] - default argument
*  @param LS - Lsystem structure containing the alphabet
*  @param individuals - array with genomes
**/
void Evolution::developIndividuals(
    int argc,
    char *argv[],
    LSystem LS,
    int generation,
    std::vector< Genome > &individuals,
    std::string dir)
{
  // for each genome in the array
  for (size_t i = 0; i < individuals.size(); ++i)
  {
    // develops genome
    individuals[i].developGenomeIndirect(
        argc,
        argv,
        this->params,
        LS,
        generation,
        this->path+"experiments/"+dir);
  }
}

/**
 * Loads population of genomes from files, from previous experiment.
 **/
void Evolution::loadIndividuals(int generation, std::string type)
{

  std::string path_list = "";
  std::string folder = "";

  if(generation==1)
    folder = "offspringpop";
  else
    folder = "selectedpop";

  if(type == "population")
  {
    // generates list of files (genomes of last population)
    std::system(("ls " + this->path + "experiments/" + this->experiment_name +
                 "/"+folder + std::to_string(generation) +
                 ">" + this->path + "experiments/" + this->experiment_name +
                 "/temp.txt").c_str());

    path_list = this->path + "experiments/" +
                this->experiment_name +"/temp.txt";
  }
  else
  {
    // reads list of genomes in archive
    path_list = this->path+"experiments/" +
                this->experiment_name + "/archive.txt";
  }

  std::ifstream listgenomes(path_list.c_str());
  std::string linegenome;

  //for each file (genome)
  while (getline(
      listgenomes,
      linegenome))
  {

    std::vector< std::string > tokens;
    boost::split(
        tokens,
        linegenome,
        boost::is_any_of("_"));

    // if it is a body file for pop recovery, or archive recovery
    if(tokens[0] == "body" or type != "population")
    {
      std::string idgenome = "";
      std::string idparent1 = "";
      std::string idparent2 = "";

      if (type == "population")
      {
        boost::split(
            tokens,
            linegenome,
            boost::is_any_of("_."));
        idgenome = tokens[1];
        idparent1 = tokens[3];
        idparent2 = tokens[5];
      }
      else
      {
        boost::split(
            tokens,
            linegenome,
            boost::is_any_of(" "));
        idgenome = tokens[0];
        idparent1 = tokens[1];
        idparent2 = tokens[2];
      }

      Genome gen = Genome(
          idgenome,
          idparent1,
          idparent2);

      // finds number of generation to which the genome belongs to
      int generation_genome = this->getGeneration_genome(idgenome);

      // reads the file with the genome
      std::ifstream listalphabet(
          this->path + "experiments/" + this->experiment_name +
          "/offspringpop" +
          std::to_string(generation_genome) + "/genome" + idgenome +
          ".txt");
      std::string linealphabet;
      // for each letter of the alphabet
      while (getline(
          listalphabet,
          linealphabet))
      {

        // gets letter and production rule from file
        std::vector< std::string > items;
        boost::split(
            items,
            linealphabet,
            boost::is_any_of(" "));
        std::vector< std::string > items_rule(
            items.begin() + 1,
            items.begin() + items.size() -
            1);

        // build a genetic-string with the production rule for the letter
        auto lgs = GeneticString();
        lgs = gen.build_genetic_string(
            lgs,
            items_rule);

        // adds letter and its production rule (made a genetic-string) to the grammar of the genome
        gen.addLetterGrammar(
            items[0],
            lgs);

      }

      // reads the measures of the genome
      std::ifstream listmeasures(
          this->path + "experiments/" + this->experiment_name +
          "/offspringpop" +
          std::to_string(generation_genome) + "/measures" + idgenome +
          ".txt");
      std::string linemeasures;
      // for each measure of the list
      while (getline(
          listmeasures,
          linemeasures))
      {

        std::vector< std::string > tokens;
        boost::split(
            tokens,
            linemeasures,
            boost::is_any_of(":"));

        gen.updateMeasure(
            tokens[0],
            std::stod(tokens[1]));
      }

      // reads fitness of the genome
      std::ifstream fitness(
          this->path + "experiments/" + this->experiment_name +
          "/offspringpop" +
          std::to_string(generation_genome) + "/fitness_" + idgenome +
          ".txt");
      if (fitness.is_open())
      {
        std::string linefitness;
        getline(
            fitness,
            linefitness);
        gen.updateLocomotionFitness(std::stod(linefitness));
      }

      std::ifstream fitness2(
          this->path + "experiments/" + this->experiment_name +
          "/offspringpop" +
          std::to_string(generation_genome) + "/fitness2_" + idgenome +
          ".txt");
      if (fitness2.is_open())
      {
        std::string linefitness;
        getline(
            fitness2,
            linefitness);
        gen.updateBalanceFitness(std::stod(linefitness));
      }


      if (type == "population")
      {
        // adds genome to the population
        this->population.push_back(gen);
      }
      else
      {
        // adds genome to the archive
        this->archive.push_back(gen);
      }
    }

  }

};


/**
 *  Loads state of previous experiment.
 **/
int Evolution::loadExperiment()
{
  // loads state of parameters from previous experiment
  this->loadsParams();

  this->aux = Aux(
      this->experiment_name,
      this->getParams(),
      this->path);

  this->logsTime("start recovery gen");

  // loads generation number from previous  experiment
  int gi = std::stoi(this->readsEvolutionState()[0]);
  // loads next_id from previous experiment
  this->next_id = std::stoi(this->readsEvolutionState()[1]);

  // deletes possible remains of unfinished generation
  std::string pathdir =
      this->path+"experiments/" + this->experiment_name + "/selectedpop" +
      std::to_string(gi + 1);
  system(("exec rm -r " + pathdir).c_str());
  pathdir = this->path+"experiments/" + this->experiment_name + "/offspringpop" +
            std::to_string(gi + 1);
  system(("exec rm -r " + pathdir).c_str());

  std::cout<<"clean folders "<<std::endl;
  // loads  population and archive
  this->loadIndividuals(gi, "population");

  std::cout<<"load pop "<<std::endl;
  this->loadIndividuals(0, "archive");

  std::cout<<"clean archive "<<std::endl;
  // loads state of the morphological_grid_accumulated
  std::string line;
  std::ifstream myfile(
      this->path+"experiments/" + this->experiment_name +
      "/morphological_grid_accumulated.txt");
  while (getline(
      myfile,
      line))
  {
    std::vector< std::string > tokens, tokens2;
    boost::split(
        tokens,
        line,
        boost::is_any_of("-"));
    boost::split(
        tokens2,
        tokens[1],
        boost::is_any_of(" "));
    std::vector< std::string > tokens3(
        tokens2.begin(),
        tokens2.begin() + tokens2.size() -
        1);
    std::vector< std::string > points;
    for (int i = 0;
         i < tokens3.size();
         i++)
    {
      points.push_back(tokens3[i]);
    }
    this->morphological_grid_accumulated[tokens[0]] = points;
  }
  std::cout<<"load grid "<<std::endl;
  myfile.close();

  return gi;
}

/* Saves the fitness for genome.
 * */
void Evolution::saveLocomotionFitness(
    std::string genome_id,
    double fitness)
{

  // updates locomotion fitness for the genome
  int index=0;
  while(this->offspring[index].getId() != genome_id)
    index++;
  this->offspring[index].updateLocomotionFitness(fitness);

  int generation_genome = this->getGeneration_genome(this->offspring[index].getId());

  // exports fitness for recovery
  std::ofstream file;
  std::string path2 =
      this->path+"experiments/"
      + this->experiment_name + "/offspringpop" +
      std::to_string(generation_genome)+ "/fitness_"+this->offspring[index]
          .getId()+".txt";
  file.open(path2);
  file
      <<  this->offspring[index].getLocomotionFitness();
  file.close();
}

/* Saves the fitness for genome.
 * */
void Evolution::saveBalanceFitness(
    std::string genome_id,
    double fitness)
{

  // updates balance fitness for the genome
  int index=0;
  while(this->offspring[index].getId() != genome_id)
    index++;
  this->offspring[index].updateBalanceFitness(fitness);

  int generation_genome = this->getGeneration_genome(this->offspring[index].getId());

  // exports fitness for recovery
  std::ofstream file;
  std::string path2 =
      this->path+"experiments/"
      + this->experiment_name + "/offspringpop" +
      std::to_string(generation_genome)+ "/fitness2_"+this->offspring[index]
          .getId()+".txt";
  file.open(path2);
  file
      <<  this->offspring[index].getBalanceFitness();
  file.close();
}


/**
 * Tries to add individuals to an archive for NS.
 **/
void Evolution::addToArchive()
{
  std::random_device rd;
  // distribution for 0-1 probabilities
  std::default_random_engine generator(rd());
  std::uniform_real_distribution< double > prob(0.0, 1.0);

  std::ofstream file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/archive.txt";
  file.open(path, std::ofstream::app);

  for (int i = 0; i < this->offspring.size(); i++)
  {
    // if raffled probability is within the constrained probability
    if (prob(generator) < this->getParams()["prob_add_archive"])
    {
      //copies object of the genome to archive
      this->archive.push_back(this->offspring[i]);

      file <<  this->offspring[i].getId()
           << " " << this->offspring[i].getId_parent1()
           << " " << this->offspring[i].getId_parent2()<< std::endl;
    }
  }

  file.close();
}



/* Calculate the novelty of the individuals.
 * */
void Evolution::calculateNovelty()
{
  std::vector<Genome> individuals_compare;
  individuals_compare.insert(individuals_compare.end(),
                             this->population.begin(), this->population.end());
  individuals_compare.insert(individuals_compare.end(),
                             this->archive.begin(), this->archive.end());

  //matrix with all individuals
  // columns: number of metrics / lines: number of genomes
  arma::mat compare(
      individuals_compare[0].getMeasures().size(),
      individuals_compare.size());


  for (int i = 0; i < individuals_compare.size(); i++)
  {
    int m = 0;
    for (const auto &it : individuals_compare[i].getMeasures())
    {
      compare(
          m,
          i) = it.second;
      m++;
    }
  }

  for (int i = 0; i < this->population.size(); i++)
  {
    // matrix with individuals which will be compared to the others
    // columns: number of metrics / single line: genome
    arma::mat reference(
        this->population[0].getMeasures().size(),
        1);

    int m = 0;
    for (const auto &it : this->population[i].getMeasures())
    {

      reference(
          m,
          0) = it.second;
      m++;
    }

    NeighborSearch< NearestNeighborSort, EuclideanDistance > nn(compare);
    arma::Mat< size_t > neighbors;
    arma::mat distances;

    // search for each individual, the nearest neighbors (+1 because it includes itself)
    nn.Search(
        reference,
        this->params["k_neighbors"] + 1,
        neighbors,
        distances);

    double fitness = 0;
    for (size_t j = 0; j < neighbors.n_elem; ++j)
    {
      fitness += distances[j];
      this->aux.logs(
          "nearest neighbor  " + std::to_string(j) + " for genome "
          + this->population[i].getId() + " has distance "
          + std::to_string(distances[j]));
    }

    // averages the nearest neighboards
    fitness = fitness / this->params["k_neighbors"];

    this->population[i].updateNoveltyFitness(fitness);

  }
}

void Evolution::calculateNoveltyLocomotion()
{
  std::vector<Genome> individuals_compare;
  individuals_compare.insert(individuals_compare.end(),
                             this->population.begin(), this->population.end());
  individuals_compare.insert(individuals_compare.end(),
                             this->archive.begin(), this->archive.end());

  //matrix with all individuals
  // columns: number of metrics / lines: number of genomes
  arma::mat compare(
      1,
      individuals_compare.size());

  for (int i = 0; i < individuals_compare.size(); i++)
  {
    int m = 0;
    compare(
          m,
          i) = individuals_compare[i].getLocomotionFitness();
  }

  for (int i = 0; i < this->population.size(); i++)
  {
    // matrix with individuals which will be compared to the others
    // columns: number of metrics / single line: genome
    arma::mat reference(
        1,
        1);

    int m = 0;
    reference(
          m,
          0) = this->population[i].getLocomotionFitness();


    NeighborSearch< NearestNeighborSort, EuclideanDistance > nn(compare);
    arma::Mat< size_t > neighbors;
    arma::mat distances;

    // search for each individual, the nearest neighbors (+1 because it includes itself)
    nn.Search(
        reference,
        this->params["k_neighbors"] + 1,
        neighbors,
        distances);

    double fitness = 0;
    for (size_t j = 0; j < neighbors.n_elem; ++j)
    {
      fitness += log10(distances[j])*-1;
      this->aux.logs(
          "nearest neighbor  " + std::to_string(j) + " for genome "
          + this->population[i].getId() + " has distance "
          + std::to_string(distances[j]));
    }

    // averages the nearest neighboards
    fitness = fitness / this->params["k_neighbors"];

    this->population[i].updateNoveltyLocomotionFitness(fitness);

  }
}

/**
 * Consolidates final fitness for evolution.
 **/
void Evolution::calculateFinalFitness()
{
  for (int i = 0; i < this->population.size(); i++)
  {

    double fitness = //this->population[i].getLocomotionFitness()
                   //*
        this->population[i].getNoveltyFitness()
    ;

    this->population[i].updateFinalFitness(fitness);
  }
}

void Evolution::calculatePenaltyFitness()
{
  for (int i = 0; i < this->population.size(); i++)
  {

    double fitness =
        std::max(0.1, 1 - this->population[i].getMeasures()["connectivity2"]);

    this->population[i].updatePenaltyFitness(fitness);
  }
}

/**
*  Evolution in the search for locomotion - part 1 of the process.
**/
double Evolution::runExperiment_part1(
    int generation, int load_experiment)
{
  int argc = 1;
  char *argv[] = {"a"};


  if(load_experiment == 1)
  {
    this->loadExperiment();
  }


  this->aux = Aux(
      this->experiment_name,
      this->getParams(),
      this->path);

  // loads alphabet with letters and commands
  LSystem LS(this->getParams());

  this->aux.logs("------------ generation " + std::to_string(generation) + " ------------");
  this->logsTime("start gen");

  this->aux.createFolder(this->path+"experiments/"+this->experiment_name +
                         "/offspringpop"+std::to_string(generation));

  if(generation == 1)
  {
    this->createHeader();

    // initializes population
    this->initPopulation(LS);
  }
  else{

    this->aux.createFolder(
        this->path+"experiments/"+this->experiment_name + "/selectedpop" + std::to_string(generation));

    this->offspring = std::vector< Genome >();

    // creates offspring
    this->crossover(LS);

    // mutates new individuals
    this->mutation(LS);

  }

  // develops genomes of the initial population
  this->developIndividuals(
      argc,
      argv,
      LS,
      generation,
      this->offspring,
      this->experiment_name + "/offspringpop");

  // measures phenotypes of the individuals
  this->measureIndividuals(
      generation,
      this->offspring,
      "/offspringpop");

  // updates the average measures for the population
  this->savesValidity(generation);

}


/**
*  Evolution in the search for locomotion - part 2 of the process. After
 *  fitness evaluation.
**/
double Evolution::runExperiment_part2(int generation)
{

  // adds new individuals to population
  for (int j = 0;
       j < this->offspring.size();
       j++)
  {
    this->population.push_back(this->offspring[j]);
  }

  this->calculateNovelty();

  //this->calculateNoveltyLocomotion();

  //this->calculatePenaltyFitness();

  this->calculateFinalFitness();

  // saves a history with the fitnesses of new individuals
  this->saveHistory(generation);

  std::vector< int > niche_measures = this->calculateNicheCoverage();

  if(generation != 1)
  {
    // selects individuals, keeping the population with a fixed size
    this->selection();

    // saves phenotypes of the selected population to a separated folder (only for visualization issues)
    this->exportPop(generation);
  }

  // saves metrics of evolution to file
  this->exportGenerationMetrics(
      generation,niche_measures);

  this->summaryNicheCoverage();

  // adds some of the new individuals to archive
  // fitness in the archive might be outdated by it would be used
  this->addToArchive();

  // saves the number of the last generation created/evaluated
  this->writesEvolutionState(
      generation,
      this->next_id);

  this->logsTime("end gen");
}

void Evolution::summaryNicheCoverage()
{

  std::ofstream file;

  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/morphological_grid_summary.txt";
  file.open(path);
  file << "point count"<< std::endl;

  for (const auto &it : this->morphological_grid_accumulated)
  {

    file << it.first + " " << it.second.size() << std::endl;
  }
  file.close();

}


std::vector< Genome > Evolution::getPopulation()
{
  return this->population;
}

std::map< std::string, double > Evolution::getParams()
{
  return params;
}


/*
 * Exports the parameters of the experiment.
 **/
void Evolution::saveParameters()
{

  std::ofstream param_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/configuration.txt";
  param_file.open(path);

  // writes each parameter to a different line in a the file
  for (auto &it : this->getParams())
  {

    param_file << it.first << " " << it.second;
    param_file << std::endl;
  }
  param_file.close();
}


/*
 * Logs time.
 **/
void Evolution::logsTime(std::string moment)
{

  time_t sta = time(0);
  char *dtsta = ctime(&sta);
  this->aux.logs("experiment " + moment + ": " + dtsta);

}

/*
 * Logs a reference of generation and genome for evolution to be recovered.
 * */
void Evolution::writesEvolutionState(
    int generation,
    int next_id)
{

  std::ofstream logs_file;
  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/evolutionstate.txt";
  logs_file.open(path);
  logs_file << generation << " " << next_id;
  logs_file.close();
}

/*
 * Reads number of the generation from which the recovered evolution should start from.
 * */
std::vector< std::string > Evolution::readsEvolutionState()
{

  std::string line;
  std::ifstream myfile(
      this->path+"experiments/" + this->experiment_name +
      "/evolutionstate.txt");
  if (myfile.is_open())
  {

    getline(
        myfile,
        line);
    std::vector< std::string > tokens;
    // parameters label and value separated by space
    boost::split(
        tokens,
        line,
        boost::is_any_of(" "));

    return tokens;

  }
  else
  {
    this->aux.logs("Unable to open evolutionstate file.");
  }
  myfile.close();

}


/*
 * Calculates the quality metric for the novelty search: niche coverage and bins of measures
 * */
std::vector< int >
Evolution::calculateNicheCoverage()
{
  morphological_grid_generation =  std::map<std::string, std::vector<double>>();

  for (int i = 0;
       i < this->offspring.size();
       i++)
  {

    std::string key_point = "";
    double distance = 0;

    // for each measure (dimension)
    for (const auto &it : this->offspring[i].getMeasures())
    {
      // accounts for NC
      // for each bin
      for (int b = 1;
           b <= this->params["grid_bins"];
           b++)
      {
        // if value is zero, sets into the first bin
        if (it.second == 0 and b == 1)
        {

          key_point += std::to_string(b) + "|";
          distance +=
              -1 * (it.second - (b / this->params["grid_bins"]));
        }
        // otherwise, sets value for measure into the correct bin
        if (it.second > ((b - 1) / this->params["grid_bins"]) and
            it.second <= (b / this->params["grid_bins"]))
        {
          key_point += std::to_string(b) + "|";
          distance +=
              -1 * (it.second - (b / this->params["grid_bins"]));
        }
      }
    }

    // if point already exists in the array, adds an individual and the difference between them
    if(morphological_grid_generation.count(key_point)>0) {

      morphological_grid_generation[key_point].push_back(distance); // add map with key=id value=distance ?

      // if point does not exist in the array yet, , adds new point with its first individual and the difference between them
    }else {
      std::vector<double> individual; individual.push_back(distance);
      morphological_grid_generation[key_point] = individual;
    }

    // if point already exists in the array, adds an individual and the id
    if (this->morphological_grid_accumulated.count(key_point) > 0)
    {

      this->morphological_grid_accumulated[key_point].push_back
          (this->offspring[i].getId());

      // if point does not exist in the array yet, , adds new point with its first individual and the difference between them
    }
    else
    {
      std::vector< std::string > individual;
      individual.push_back(this->offspring[i].getId());
      this->morphological_grid_accumulated[key_point] = individual;
    }
  }



  // logs state of the grid
  std::ofstream myfile;
  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/morphological_grid_accumulated.txt";
  myfile.open(path);
  for (const auto &it : this->morphological_grid_accumulated)
  {
    myfile << it.first + "-";
    for (int i = 0;
         i < it.second.size();
         i++)
    {
      myfile << it.second[i] << " ";
    }
    myfile << std::endl;
  }
  myfile.close();


  std::vector< int > morphological_grids;
  morphological_grids.push_back((int)morphological_grid_generation.size());
  morphological_grids.push_back((int) this->morphological_grid_accumulated.size());


  return morphological_grids;
}