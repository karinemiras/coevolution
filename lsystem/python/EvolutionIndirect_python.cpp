
#include <boost/python.hpp>
#include "EvolutionIndirect.h"

class EvolutionIndirect_python : public EvolutionIndirect {
public:
  EvolutionIndirect_python(std::string experiment_name,
                      std::string type_experiment)
     : EvolutionIndirect(experiment_name, type_experiment)
  {}


//  double runExperiment()
//  {
//    EvolutionIndirect::runExperiment(0, nullptr);
//  }

};

BOOST_PYTHON_MODULE (lsystem_python)
{
  boost::python::class_< EvolutionIndirect_python, boost::noncopyable >(
          "EvolutionIndirect",
          boost::python::init< std::string,
                               std::string >())
      .def("setupEvolution",
           &EvolutionIndirect_python::setupEvolution)
      .def("runExperiment_part1",
           &EvolutionIndirect_python::runExperiment_part1)
      .def("runExperiment_part2",
           &EvolutionIndirect_python::runExperiment_part2)
      .def("saveLocomotionFitness",
           &EvolutionIndirect_python::saveLocomotionFitness)
      .def("saveBalanceFitness",
           &EvolutionIndirect_python::saveBalanceFitness)
      ;
}
