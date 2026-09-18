#ifndef PARAMS_H
#define PARAMS_H

#include "circuit.h"

namespace snrk {

inline int SplinePartition = 3;

inline in::witness_t wStart = 1;
inline in::witness_t wStep = 1;

inline bool MultiThreading = true;

#define GateFuncHeader (const snrk::CircuitValue &a, const snrk::CircuitValue &b, const snrk::CircuitValue &c)
using GateFunc = std::function<std::vector<BaseGate>GateFuncHeader>;
#define GateLambdaHeader GateFuncHeader -> std::vector<snrk::BaseGate>

inline std::unordered_map<snrk::GateType, GateFunc> Gates;

}
#endif // PARAMS_H
