#ifndef PROOF_H
#define PROOF_H

#include "setup.h"

namespace snrk {

class ProverProof : public Jsonable
{
public:
    ProverProof() = default;

    ProverProof(const CircutParams &cp, ProverParams &pp, const std::string &trInit);

    virtual bool check(const VerifierParams &vp, const std::string &trInit);

public:
    virtual json_t toJson() const override;

    virtual bool fromJson(const json_t &json) override;

private:
    in::SplinePolynom correctGates(const CircutParams::TParams_t &TParams, const CircutParams::SParams_t &SParams) const;

    in::SplinePolynom correctPermulations(const in::witnesses_t &witnesses, const CircutParams::TParams_t &TParams,
                                          const CircutParams::WParams_t &WParams, std::array<in::value_t, 2> trPow) const;

    MerkleTree makeTree(const in::SplinePolynom &polynom, ProverParams &pp, in::Transcript &tr) const;

    in::commit_t m_Q, m_F;

    in::value_t m_Pr, m_QPr, m_Gr, m_QGr;
};

}

#endif // PROOF_H
