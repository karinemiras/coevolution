//
// Created by Karine Miras on 15/05/2017.
//

#ifndef LSYSTEM_PROTO_AUX_H
#define LSYSTEM_PROTO_AUX_H

#include <map>
#include <string>

class Aux{


    public:

    Aux(std::string experiment_name,
        std::map<std::string, double> params,
        std::string path)
    {
        this->params = params;
        this->experiment_name = experiment_name;
        this->path = path;
    }

    void createFolder(std::string dir);
    void removeFolder(std::string dir);
    void logs( std::string message);

    std::string experiment_name;
    std::map<std::string, double> params;
    std::string path;
};


#endif //LSYSTEM_PROTO_AUX_H
