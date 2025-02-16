#include "Alpha.H"

using namespace Foam;

preciceAdapter::FF::Alpha::Alpha(
        const Foam::fvMesh &mesh,
        const std::string nameAlpha,
        DimensionReduction &red)
        : Alpha_(
        const_cast<volScalarField *>(
                &mesh.lookupObject<volScalarField>(nameAlpha))),
          red{red} {
    dataType_ = scalar;
}

std::size_t preciceAdapter::FF::Alpha::write(double *buffer, bool meshConnectivity, const unsigned int dim) {
    int bufferIndex = 0;


    if (this->locationType_ == LocationType::volumeCenters) {
        if (cellSetNames_.empty()) {
            for (const auto &cell: Alpha_->internalField()) {
                buffer[bufferIndex++] = cell;
            }
        } else {
            for (const auto &cellSetName: cellSetNames_) {
                cellSet overlapRegion(Alpha_->mesh(), cellSetName);
                const labelList &cells = overlapRegion.toc();

                for (const auto &currentCell: cells) {
                    // Copy the alpha valus into the buffer
                    buffer[bufferIndex++] = Alpha_->internalField()[currentCell];
                }
            }
        }
    }

    auto heights = red.getH();
    std::cout << "Heights " << heights.size() << "\n";
    for (auto h: heights) {
        buffer[bufferIndex++] = h;
    }

    return bufferIndex;
}

void preciceAdapter::FF::Alpha::read(double *buffer, const unsigned int dim) {
    int bufferIndex = 0;

    if (this->locationType_ == LocationType::volumeCenters) {
        if (cellSetNames_.empty()) {
            for (auto &cell: Alpha_->ref()) {
                cell = buffer[bufferIndex++];
            }
        } else {
            for (const auto &cellSetName: cellSetNames_) {
                cellSet overlapRegion(Alpha_->mesh(), cellSetName);
                const labelList &cells = overlapRegion.toc();

                for (const auto &currentCell: cells) {
                    // Copy the pressure into the buffer
                    Alpha_->ref()[currentCell] = buffer[bufferIndex++];
                }
            }
        }
    }

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++) {
        int patchID = patchIDs_.at(j);
        auto &mesh = Alpha_->mesh();
        auto &patch_cells = Alpha_->mesh().boundary()[patchID].faceCells();
        auto &cells = mesh.cells();

        // For every water column get its height
        int z = 0;
        for (const auto &c_label: patch_cells) {
            auto height = buffer[bufferIndex++];
            auto c = cells[c_label];

            // Calculate water height
            auto bounding_box = c.box(mesh);
            double top = max(bounding_box.first().z(), bounding_box.second().z());
            double bottom = min(bounding_box.first().z(), bounding_box.second().z());

//            if (top < height) {
//                Alpha_->boundaryFieldRef()[patchID][z] = 1;
//            } else if (bottom >= height) {
//                Alpha_->boundaryFieldRef()[patchID][z] = 0;
//            } else {
//                Alpha_->boundaryFieldRef()[patchID][z] = (height - bottom) / (top - bottom);
//            }

            if (top < height) {
                Alpha_->ref()[c_label] = 1;
            } else if (bottom >= height) {
                Alpha_->ref()[c_label] = 0;
            } else {
                Alpha_->ref()[c_label] = (height - bottom) / (top - bottom);
            }

            z++;
        }
    }
}

bool preciceAdapter::FF::Alpha::isLocationTypeSupported(const bool meshConnectivity) const {
    if (meshConnectivity) {
        return (this->locationType_ == LocationType::faceCenters);
    } else {
        return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::volumeCenters);
    }
}

std::string preciceAdapter::FF::Alpha::getDataName() const {
    return "Alpha";
}
