//
// Created by armin on 11/25/24.
//

#include "DimensionReduction.h"

#include <Adapter.H>

#include "uniformSet.H"
#include "meshSearch.H"

preciceAdapter::FF::DimensionReduction::DimensionReduction(const Foam::fvMesh &mesh,
                                                           const std::string nameAlpha,
                                                           const std::string nameV,
                                                           std::vector<int> patchIDs,
                                                           int n_samples) :
        Alpha(const_cast<volScalarField *>(&mesh.lookupObject<volScalarField>(nameAlpha))),
        V(const_cast<volVectorField *>(&mesh.lookupObject<volVectorField>(nameV))),
        n_samples{n_samples},
        patchIDs{patchIDs},
        nameAlpha{nameAlpha} {
    meshSearch searchEngine(mesh);

    adapterInfo("Start Dimension Mapping");
    int counter = 0;

    // For every boundary patch of the interface
    for (int patchID: patchIDs) {
        auto &temp_mesh = Alpha->mesh();
        auto &face_centers = temp_mesh.boundary()[patchID].Cf();

        // For every water column get its height and average water velocity
        for (const auto &center: face_centers) {

            sample_size = center.z() / n_samples;
            double offset = sample_size / 2;
            point start(center.x(), center.y(), center.z() - offset);
            point end(center.x(), center.y(), offset);

//            Foam::dictionary setdict{};
//            setdict.add("type", "uniform");
//            setdict.add("axis", "z");
//            setdict.add("start", "(" + std::to_string(start.x()) + " " + std::to_string(start.y()) + " " + std::to_string(start.z()) + ")");
//            setdict.add("end", "(" + std::to_string(end.x()) + " " + std::to_string(end.y()) + " " + std::to_string(end.z()) + ")");
//            setdict.add("nPoints", std::to_string(n_samples));

            sampleSets.emplace_back("Sample" + std::to_string(counter++), temp_mesh, searchEngine, "z", start, end,
                                    n_samples);

//            for (auto cell: sampleSets[sampleSets.size() - 1].cells()) {
//                auto p = mesh.cellCentres()[cell];
//                std::cout << "(" << p.x() << "," << p.y() << "," << p.z() << ")\n";
//            }
        }
    }

    adapterInfo("Finish Dimension Mapping");
}

std::vector<double> &preciceAdapter::FF::DimensionReduction::getH() {
    if (!calculated || read_H) calculateFields();
    read_H = true;
    if (read_H && read_V) calculated = false;

    return heights;
}

void preciceAdapter::FF::DimensionReduction::calculateFields() {
    DEBUG(adapterInfo("Start 3D-2D reduction calculation"));

    const double WATER_LEVEL = 1000;

    velocities.clear();
    heights.clear();

    for (uniformSet &set: sampleSets) {
        double alpha_sum = 0;
        double u_sum = 0;
        double v_sum = 0;

        bool obstacle = false;

        for (int i = 0; i < set.cells().size(); i++) {
            auto cell = set.cells()[i];
            alpha_sum += Alpha->internalField()[cell];
            u_sum += V->internalField()[cell].x() * Alpha->internalField()[cell];
            v_sum += V->internalField()[cell].y() * Alpha->internalField()[cell];
        }

        double calc_height = alpha_sum * sample_size;

//        if (obstacle && calc_height < WATER_LEVEL)
        if (set.cells().size() < n_samples) {
            calc_height = WATER_LEVEL;
        }

        heights.push_back(calc_height);
        if (u_sum == 0) velocities.push_back(0);
        else velocities.push_back(u_sum / alpha_sum);
        if (v_sum == 0) velocities.push_back(0);
        else velocities.push_back(v_sum / alpha_sum);
        velocities.push_back(0); // Precice expects 3D vector
    }

    calculated = true;
    read_H = false;
    read_V = false;

    DEBUG(adapterInfo("Finish 3D-2D reduction"));
}

std::vector<double> &preciceAdapter::FF::DimensionReduction::getV() {
    if (!calculated || read_V) calculateFields();
    read_V = true;
    if (read_H && read_V) calculated = false;

    return velocities;
}
