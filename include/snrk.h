#ifndef SNARK_H
#define SNARK_H

#include "circuit.h"
#include "proof.h"

namespace snrk {

inline void init()
{
    mcl::initPairing(mcl::BN254);

    MultiThreading = true;

    SplinePartition = 6;
}

}

#endif // SNARK_H
