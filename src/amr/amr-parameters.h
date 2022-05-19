#pragma once

typedef struct amrParameters {
    unsigned refinementFactor;
    unsigned refinementInterval;
    unsigned bufferSize;
    double efficiencyThreshold;
    double truncationErrorTreshold;
} AMRParameters;
