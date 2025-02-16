//
// Created by armin on 10/9/24.
//

#ifndef SWEET_PRECICE_COUPLING_COUPLINGUTILITY_HPP
#define SWEET_PRECICE_COUPLING_COUPLINGUTILITY_HPP

#include <src/include/sweet/Shacks/Base.hpp>
#include <src/include/sweet/Tools/ProgramArguments.hpp>

class Coupling_Shack : public sweet::Shacks::Base {
public:
    bool use_precice = false;
    bool couple_alpha = false;

    void printProgramArguments(const std::string &i_prefix) override {
        std::cout << i_prefix << "-p ...: Use preCICE coupling\n";
    }

    bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override {
        if (i_pa.argumentWithKeyExists("-p")) {
            std::string test;
            i_pa.getArgumentValueByKey("-p", test);

            if (test == "both") {
                std::cout << "Using coupling mode BOTH\n";
                use_precice = true;
                couple_alpha = true;
            } else if (test == "simple") {
                std::cout << "Using coupling mode SIMPLE\n";
                use_precice = true;
            } else {
                throw new std::invalid_argument(test + " is not a allowed precice mode");
            }
        } else {
            std::cout << "preCICE coupling DISABLED\n";
        }

        ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa)
    }

    void printShack(const std::string &i_prefix) override {
        std::cout << "Precice: " << use_precice;
    }
};

#endif //SWEET_PRECICE_COUPLING_COUPLINGUTILITY_HPP
