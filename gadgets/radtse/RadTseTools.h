#pragma once

#include "mri_core_data.h"

namespace Gadgetron { 
    namespace RadTse {
        void newViewSelect(std::vector<std::vector<int>> &view, int nviews, int etl);

        long prepareData(double *mask, int ncol, int nviews, int etl, int echoNum, double *denscomp, double USFactor);
    }
}