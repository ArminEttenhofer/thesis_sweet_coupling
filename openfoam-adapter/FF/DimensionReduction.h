//
// Created by armin on 11/25/24.
//

#ifndef SWEET_PRECICE_COUPLING_DIMENSIONREDUCTION_H
#define SWEET_PRECICE_COUPLING_DIMENSIONREDUCTION_H

#include <meshSearch.H>
#include <uniformSet.H>
#include <vector>
#include "fvCFD.H"

namespace preciceAdapter {
    namespace FF {

        class DimensionReduction {
            bool calculated{false};
            bool read_V{false};
            bool read_H{false};

            std::vector<double> velocities{};
            std::vector<double> heights{};

            Foam::volScalarField *Alpha{};
            Foam::volVectorField *V{};
            std::vector<Foam::uniformSet> sampleSets{};
            int n_samples;
            double sample_size{};

            std::vector<int> patchIDs;

            void calculateFields();

        public:
            std::string nameAlpha{};

            DimensionReduction(const Foam::fvMesh &mesh,
                               const std::string nameAlpha,
                               const std::string nameV,
                               std::vector<int> patchIDs,
                               int n_samples);

            std::vector<double> &getV();

            std::vector<double> &getH();
        };

    }
}


#endif //SWEET_PRECICE_COUPLING_DIMENSIONREDUCTION_H
