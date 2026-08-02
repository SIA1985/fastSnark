#ifndef SETUP_H
#define SETUP_H

#include "polynom.h"

#include <unordered_set>

namespace in {
    struct SplittedT_t
    {
        InterpolationPolynom left;
        InterpolationPolynom right;
        InterpolationPolynom result;
    };

    witnesses_t genWitnesses(in::witness_t start, std::size_t count, in::witness_t wStep);

}

namespace snrk {

class GlobalParams
{
public:
    GlobalParams(in::value_t t);

    ProverParams PP() const;

    VerifierParams VP() const;

private:
    ProverParams m_pp;
    VerifierParams m_vp;
};

void to_json(json_t& j, const GlobalParams& gp);
void from_json(const snrk::json_t& j, GlobalParams& gp);

class CircutParams {
public:
    using opsFromS_t = std::unordered_map<BaseGateType, in::dots_t>;
    using TParams_t = struct{in::T_t t; in::SplittedT_t splittedT;};
    using SParams_t = struct{opsFromS_t opsFromS;};
    using WParams_t = struct{in::W_t wt; in::W_t wi;};

    using cond_t = std::unordered_set<in::witness_t>;

    using params_t = struct{TParams_t TParams; SParams_t SParams; WParams_t WParams;};


    CircutParams(Circuit &circuit, ProverParams &pp);

    in::witnesses_t witnesses() const;

    std::size_t witnessesCount() const;

    params_t params() const;

    std::size_t circuitSize() const;

private:
    void generateT(Circuit &circuit);

    void generateS(Circuit &circuit);

    void generateW(Circuit &circuit);

    in::witnesses_t m_witnesses;

    ProverParams m_pp;

    in::T_t m_T;
    in::SplittedT_t m_splittedT;

    opsFromS_t m_opsFromS;

    in::W_t m_WT, m_WI;

    std::size_t m_circuitSize;
};

}

#endif // SETUP_H
