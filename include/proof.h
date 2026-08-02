#ifndef PROOF_H
#define PROOF_H

#include "setup.h"

namespace in {

class Transcript {
public:
   Transcript(const std::string &init);

   void appendHash(const std::string &label, const hash_t &hashData);

   void appendPoint(const std::string &label, const G1 &point);

   void appendScalar(const std::string &label, const value_t &scalar);

   value_t operator()(const std::string &label);

   value_t getScaled(const std::string &label);

   template<std::size_t N>
   std::array<value_t, N> getPowerArray(const std::string &label);

private:
   void updateState(const std::string &data);

   hash_t m_state{""};
};

}

namespace snrk {

class ProverProof : public Jsonable
{
    using WResult_t = struct{in::SplinePolynom W, WShift1;};

    using BatchCommit_t = struct {
    in::commit_t T, WT, WI, W, WNext, WitnessZ,
                 QP, QG, QC, G, C, PartitionZ;
    };

    using BatchY_t = struct {
    in::Y_t T, WT, WI, W, WNext, WitnessZ,
            QP, QG, QC, G, C, PartitionZ;
    };

public:
    ProverProof() = default;

    ProverProof(const CircutParams &cp, ProverParams &pp, const std::string &trInit);

    virtual bool check(const VerifierParams &vp, const std::string &trInit);

public:
    virtual json_t toJson() const override;

    virtual bool fromJson(const json_t &json) override;

private:
    in::SplinePolynom correctGates(const in::SplittedT_t &t, const CircutParams::SParams_t &SParams) const;

    WResult_t correctPermulations(const in::witnesses_t &witnesses, const in::SplinePolynom &num, const in::SplinePolynom &den) const;

    in::SplinePolynom correctСontinuity(in::SplinePolynom poly) const;

    in::hash_t rWitnessLeaf(int i, std::array<in::value_t, 12> alphaPow) const;

    MerkleTree makeTree(in::SplinePolynom &polynom, ProverParams &pp) const;

    in::hash_t m_witnessMerkleRoot;
    MerkleTree::MultiProof_t m_witnessMerkleProof;

    std::vector<BatchCommit_t> m_batchCommits;
    std::vector<BatchY_t> m_batchY;

    /*Защита от полиномов высокой степени*/
    in::commits_t m_piLow, m_piMid, m_piHigh;
};

}

#endif // PROOF_H
