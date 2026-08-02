#include "proof.h"

#include "threader.h"


namespace in {

Transcript::Transcript(const std::string &init)
{
    updateState(init);
}

void Transcript::appendHash(const std::string &label, const hash_t &hashData)
{
    updateState(label + hashData);
}

void Transcript::appendPoint(const std::string &label, const G1 &point)
{
    char buf[1024];
    size_t n = point.serialize(buf, sizeof(buf));
    updateState(label + std::string{buf, n});
}

void Transcript::appendScalar(const std::string &label, const value_t &scalar)
{
    char buf[1024];
    size_t n = scalar.serialize(buf, sizeof(buf));
    updateState(label + std::string{buf, n});
}

value_t Transcript::operator()(const std::string &label)
{
    updateState(label);

    value_t challenge;
    challenge.setHashOf({m_state.begin(), m_state.end()});

    return challenge;
}

value_t Transcript::getScaled(const std::string &label)
{
    return snrk::CircuitValue::scale(this->operator()(label));
}

void Transcript::updateState(const std::string &data)
{
    m_state = hash(m_state + data);
}

template<std::size_t N>
std::array<value_t, N> Transcript::getPowerArray(const std::string &label)
{
    auto value = this->operator()(label);

    std::array<value_t, N> pow = {1};
    for(int i = 1; i < N; i++) {
        pow[i] = pow[i - 1] * value;
    }

    return pow;
}

}

namespace snrk {

//todo: auto -> тип
ProverProof::ProverProof(const CircutParams &cp, ProverParams &pp, const std::string &trInit)
{
    const auto &[TParams, SParams, WParams] = cp.params();
    const auto &witnesses = cp.witnesses();

    //Инициализация Транскрипта
    in::Transcript tr(trInit);

    //Построение всех witness-сплайнов и фиксация в дереве Меркла
    auto WitnessZ = in::ZeroWitnessPolynom(witnesses).toSplinePolynom();

    const auto &T = TParams.t.toSplinePolynom();
    const auto &WT = WParams.wt.toSplinePolynom();
    const auto &WI = WParams.wi.toSplinePolynom();

    auto alphaPow = tr.getPowerArray<12>("alpha");
    auto beta = tr("beta");
    auto gamma = tr.getScaled("gamma");

    //Проверка перестановок
    const auto &num = T + (WT * beta) + gamma;
    const auto &den = T + (WI * beta) + gamma;

    const auto &[W, WNext] = correctPermulations(witnesses, num, den);

    const auto &P = (WNext * den) - (W * num);
    const auto &QP = P / WitnessZ;

    //Проверка вычислений
    const auto &G = correctGates(TParams.splittedT, cp.params().SParams);
    const auto &QG = G / WitnessZ;

    //Проверка непрерывности сплайнов
    auto Spline = T * alphaPow[0] +
                  WT * alphaPow[1] +
                  WI * alphaPow[2] +
                  W * alphaPow[3] +
                  WNext * alphaPow[4] +
                  G * alphaPow[5] +
                  WitnessZ * alphaPow[6];

    auto PartitionZ = in::ZeroWitnessPolynom::makePartitionZeroPolynom(WitnessZ);
    const auto &C = correctСontinuity(Spline);

    const auto &QC = C / PartitionZ;
    Spline += C * alphaPow[7];

    Spline += QP * alphaPow[8];
    Spline += QG * alphaPow[9];
    Spline += QC * alphaPow[10];

    Spline += PartitionZ * alphaPow[11];
    const auto &witnessTree = makeTree(Spline, pp);

    m_witnessMerkleRoot = witnessTree.root();
    tr.appendHash("witnessMerkleRoot", m_witnessMerkleRoot);

    //Получение точек раскрытия из Тr
    in::values_t rDots(DotsForProof);
    std::vector<std::size_t> witnessIndices(DotsForProof);

    m_batchCommits.resize(DotsForProof);
    m_batchY.resize(DotsForProof);

    FORDOTS {
        auto rI = tr("r_" + std::to_string(i));
        rDots[i] = rI;

        auto rCut = in::cut(rI, witnesses.back());
        witnessIndices[i] = Spline.segment(rCut);


        const auto &rT = T[rCut];
        const auto &rWT = WT[rCut];
        const auto &rWI = WI[rCut];
        const auto &rWitnessZ = WitnessZ[rCut];
        const auto &rQG = QG[rCut];
        const auto &rW = W[rCut];
        const auto &rWNext = WNext[rCut];
        const auto &rQP = QP[rCut];
        const auto &rQC = QC[rCut];
        const auto &rG = G[rCut];
        const auto &rC = C[rCut];
        const auto &rPartitionZ = PartitionZ[rCut];

        auto &bc = m_batchCommits[i];

        bc.T = rT.commit(pp);
        bc.WT = rWT.commit(pp);
        bc.WI = rWI.commit(pp);
        bc.WitnessZ = rWitnessZ.commit(pp);
        bc.QG = rQG.commit(pp);
        bc.W = rW.commit(pp);
        bc.WNext = rWNext.commit(pp);
        bc.QP = rQP.commit(pp);
        bc.QC = rQC.commit(pp);
        bc.G = rG.commit(pp);
        bc.C = rC.commit(pp);
        bc.PartitionZ = rPartitionZ.commit(pp);

        auto &by = m_batchY[i];

        by.T = rT(rI);
        by.WT = rWT(rI);
        by.WI = rWI(rI);
        by.W = rW(rI);
        by.WNext = rWNext(rI);
        by.WitnessZ = rWitnessZ(rI);
        by.QG = rQG(rI);
        by.QP = rQP(rI);
        by.QC = rQC(rI);
        by.G = rG(rI);
        by.C = rC(rI);
        by.PartitionZ = rPartitionZ(rI);
    }

    //Получение мульти-пути Меркла
    m_witnessMerkleProof = witnessTree.multiProof(witnessIndices).value();

    //Получения доказательств
    auto delta = tr("delta");
    auto epsilon = tr("epsilon");

    m_piLow.resize(DotsForProof);
    m_piMid.resize(DotsForProof);
    m_piHigh.resize(DotsForProof);

    in::value_t currentEpsilon = 1;
    FORDOTS {
        auto rI = rDots[i];
        auto rCut = in::cut(rI, witnesses.back());
        auto ys = m_batchY[i];

        in::CanonicPolynom fLow;
        in::value_t currentDelta = 1;
        for(auto &poly : {(T[rCut] - ys.T), (WT[rCut] - ys.WT), (WI[rCut] - ys.WI),
                          (W[rCut] - ys.W), (WNext[rCut] - ys.WNext), (QC[rCut] - ys.QC),
                          (C[rCut] - ys.C), (WitnessZ[rCut] - ys.WitnessZ), (PartitionZ[rCut] - ys.PartitionZ)}) {
            fLow += poly * currentDelta;
            currentDelta *= delta;
        }

        in::CanonicPolynom fMid;
        currentDelta = 1;
        for(auto &poly : {(QG[rCut] - ys.QG), (QP[rCut] - ys.QP)}) {
            fMid += poly * currentDelta;
            currentDelta *= delta;
        }

        in::CanonicPolynom fHigh;
        currentDelta = 1;
        for(auto &poly : {(G[rCut] - ys.G)}) {
            fHigh += poly * currentDelta;
            currentDelta *= delta;
        }

        auto zI = in::CanonicPolynom({-rI, 1});


        auto qLow = (fLow / zI) * currentEpsilon;
        auto qMid = (fMid / zI) * currentEpsilon;
        auto qHigh = (fHigh / zI) * currentEpsilon;

        m_piLow[i] = qLow.commit(pp);
        m_piMid[i] = qMid.commit(pp);
        m_piHigh[i] = qHigh.commit(pp);

        currentEpsilon *= epsilon;
    }
}

bool ProverProof::check(const VerifierParams &vp, const std::string &trInit)
{
    in::Transcript tr(trInit);

    if (m_batchCommits.size() != DotsForProof ||
        m_batchY.size() != DotsForProof ||
        m_piLow.size() != DotsForProof ||
        m_piMid.size() != DotsForProof ||
        m_piHigh.size() != DotsForProof) {
        return false;
    }

    auto alphaPow = tr.getPowerArray<12>("alpha");
    auto beta = tr("beta");
    auto gamma = tr.getScaled("gamma");

    in::hashes_t witnessLeafs(DotsForProof);
    FORDOTS {
        witnessLeafs[i] = rWitnessLeaf(i, alphaPow);
    }

    if (!MerkleTree::multiVerify(m_witnessMerkleProof, witnessLeafs, m_witnessMerkleRoot)) {
        return false;
    }
    tr.appendHash("witnessMerkleRoot", m_witnessMerkleRoot);


    FORDOTS {
        const auto &by = m_batchY[i];

        in::value_t numR = by.T + (beta * by.WT) + gamma;
        in::value_t denR = by.T + (beta * by.WI) + gamma;

        in::value_t errorW = (by.WNext * denR) - (by.W * numR);

        if (errorW != by.QP * by.WitnessZ) {
            return false;
        }

        if (by.G != by.QG * by.WitnessZ) {
            return false;
        }

        if (by.C != by.QC * by.PartitionZ) {
            return false;
        }
    }


    in::values_t rDots(DotsForProof);
    FORDOTS {
        auto iStr = std::to_string(i);
        rDots[i] = tr("r_" + iStr);
    }

    auto delta = tr("delta");
    auto epsilon = tr("epsilon");
    auto dzeta = tr("dzeta");

    in::commit_t aggregatePiLow, aggregatePiLowR, aggregateFLow;
    in::commit_t aggregatePiMid, aggregatePiMidR, aggregateFMid;
    in::commit_t aggregatePiHigh, aggregatePiHighR, aggregateFHigh;

    for(auto ptr : {&aggregatePiLow, &aggregatePiLowR, &aggregateFLow,
                    &aggregatePiMid, &aggregatePiMidR, &aggregateFMid,
                    &aggregatePiHigh, &aggregatePiHighR, &aggregateFHigh}) {
        ptr->clear();
    }

    in::value_t currentEpsilon = 1;
    in::value_t currentDzeta = 1;

    FORDOTS {
        const auto &bc = m_batchCommits[i];
        const auto &by = m_batchY[i];
        auto rI = rDots[i];

        in::commit_t currentFLow;
        currentFLow.clear();
        in::value_t currentDelta = 1;
        for(auto &comm : {(bc.T - vp.g1 * by.T), (bc.WT - vp.g1 * by.WT), (bc.WI - vp.g1 * by.WI),
                          (bc.W - vp.g1 * by.W), (bc.WNext - vp.g1 * by.WNext), (bc.QC - vp.g1 * by.QC),
                          (bc.C - vp.g1 * by.C), (bc.WitnessZ - vp.g1 * by.WitnessZ), (bc.PartitionZ - vp.g1 * by.PartitionZ)}) {
            currentFLow += comm * currentDelta;
            currentDelta *= delta;
        }

        in::commit_t currentFMid;
        currentFMid.clear();
        currentDelta = 1;
        for(auto &comm : {(bc.QG - vp.g1 * by.QG), (bc.QP - vp.g1 * by.QP)}) {
            currentFMid += comm * currentDelta;
            currentDelta *= delta;
        }

        in::commit_t currentFHigh;
        currentFHigh.clear();
        currentDelta = 1;
        for(auto &comm : {(bc.G - vp.g1 * by.G)}) {
            currentFHigh += comm * currentDelta;
            currentDelta *= delta;
        }

        auto combinedWeight = currentEpsilon * currentDzeta;

        aggregateFLow += currentFLow * combinedWeight;
        aggregatePiLow += m_piLow[i] * currentDzeta;
        aggregatePiLowR += m_piLow[i] * (rI * currentDzeta);

        aggregateFMid += currentFMid * combinedWeight;
        aggregatePiMid += m_piMid[i] * currentDzeta;
        aggregatePiMidR += m_piMid[i] * (rI * currentDzeta);

        aggregateFHigh += currentFHigh * combinedWeight;
        aggregatePiHigh += m_piHigh[i] * currentDzeta;
        aggregatePiHighR += m_piHigh[i] * (rI * currentDzeta);

        currentEpsilon *= epsilon;
        currentDzeta *= dzeta;
    }

    in::GT e1Low, e2Low;
    mcl::pairing(e1Low, aggregatePiLow, vp.tG2);
    mcl::pairing(e2Low, aggregateFLow + aggregatePiLowR, vp.g2);
    if (e1Low != e2Low) {
        return false;
    }

    in::GT e1Mid, e2Mid;
    mcl::pairing(e1Mid, aggregatePiMid, vp.tG2);
    mcl::pairing(e2Mid, aggregateFMid + aggregatePiMidR, vp.g2);
    if (e1Mid != e2Mid) {
        return false;
    }

    in::GT e1High, e2High;
    mcl::pairing(e1High, aggregatePiHigh, vp.tG2);
    mcl::pairing(e2High, aggregateFHigh + aggregatePiHighR, vp.g2);
    if (e1High != e2High) {
        return false;
    }

    return true;
}

const std::string WitnessMerkleRoot = "WitnessMerkleRoot",
                  WitnessMerkleMultiProof = "WitnessMerkleMultiProof",
                  BatchCommits = "BatchCommits",
                  BatchY = "BatchY",
                  PiLow = "PiLow",
                  PiMid = "PiMid",
                  PiHigh = "PiHigh",

                  T = "T",
                  WT = "WT",
                  WI = "WI",
                  W = "W",
                  WNext = "WNext",
                  WitnessZ = "WitnessZ",
                  QP = "QP",
                  QG = "QG",
                  QC = "QC",
                  G = "G",
                  C = "C",
                  PartitionZ = "PartitionZ";

json_t ProverProof::toJson() const
{
    json_t json;

    auto fillCommits = [](snrk::json_t &json, const in::commits_t &commits)
    {
        snrk::json_t::array_t arr;
        for(const auto &commit : commits) {
            arr.push_back(toString(commit));
        }

        json = arr;
    };

    auto fillBatch = [](snrk::json_t &json, const auto &batch)
    {
        snrk::json_t::array_t arr;
        for(const auto &b : batch) {
            snrk::json_t obj;

            obj[T] = toString(b.T);
            obj[WT] = toString(b.WT);
            obj[WI] = toString(b.WI);
            obj[W] = toString(b.W);
            obj[WNext] = toString(b.WNext);
            obj[WitnessZ] = toString(b.WitnessZ);
            obj[QP] = toString(b.QP);
            obj[QG] = toString(b.QG);
            obj[QC] = toString(b.QC);
            obj[G] = toString(b.G);
            obj[C] = toString(b.C);
            obj[PartitionZ] = toString(b.PartitionZ);

            arr.push_back(obj);
        }

        json = arr;
    };

    json[WitnessMerkleRoot] = toString(m_witnessMerkleRoot);
    json[WitnessMerkleMultiProof] = m_witnessMerkleProof.toJson();
    fillBatch(json[BatchCommits], m_batchCommits);
    fillBatch(json[BatchY], m_batchY);
    fillCommits(json[PiLow], m_piLow);
    fillCommits(json[PiMid], m_piMid);
    fillCommits(json[PiHigh], m_piHigh);

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    if ((!json.contains(WitnessMerkleRoot) || !json[WitnessMerkleRoot].is_string()) &&
        (!json.contains(WitnessMerkleMultiProof) || !json[WitnessMerkleMultiProof].is_object()) &&
        (!json.contains(BatchCommits) || !json[BatchCommits].is_array()) &&
        (!json.contains(BatchY) || !json[BatchY].is_array()) &&
        (!json.contains(PiLow) || !json[PiLow].is_array()) &&
        (!json.contains(PiMid) || !json[PiMid].is_array()) &&
        (!json.contains(PiHigh) || !json[PiHigh].is_array())) {
        return false;
    }

    auto fillCommits = [](in::commits_t &commits, const snrk::json_t::array_t &items)
    {
        commits.clear();
        commits.reserve(items.size());

        for(const auto &item : items) {
            commits.push_back(fromString<in::commit_t>(item.get<std::string>()).value());
        }
    };

    auto fillBatch = [](auto &batch, const snrk::json_t::array_t &items)
    {
        using VecType = std::decay_t<decltype(batch)>;
        using BatchType = typename VecType::value_type;

        batch.clear();
        batch.reserve(items.size());

        for(const auto &item : items) {
            BatchType b;

            b.T = fromString<decltype(b.T)>(item[T].get<std::string>()).value();
            b.WT = fromString<decltype(b.WT)>(item[WT].get<std::string>()).value();
            b.WI = fromString<decltype(b.WI)>(item[WI].get<std::string>()).value();
            b.W = fromString<decltype(b.W)>(item[W].get<std::string>()).value();
            b.WNext = fromString<decltype(b.WNext)>(item[WNext].get<std::string>()).value();
            b.WitnessZ = fromString<decltype(b.WitnessZ)>(item[WitnessZ].get<std::string>()).value();
            b.QP = fromString<decltype(b.QP)>(item[QP].get<std::string>()).value();
            b.QG = fromString<decltype(b.QG)>(item[QG].get<std::string>()).value();
            b.QC = fromString<decltype(b.QC)>(item[QC].get<std::string>()).value();
            b.G = fromString<decltype(b.G)>(item[G].get<std::string>()).value();
            b.C = fromString<decltype(b.C)>(item[C].get<std::string>()).value();
            b.PartitionZ = fromString<decltype(b.PartitionZ)>(item[PartitionZ].get<std::string>()).value();

            batch.push_back(b);
        }
    };

    m_witnessMerkleRoot = fromString<in::hash_t>(json[WitnessMerkleRoot]).value();
    if (!m_witnessMerkleProof.fromJson(json[WitnessMerkleMultiProof])) {
        return false;
    }

    fillBatch(m_batchCommits, json[BatchCommits]);
    fillBatch(m_batchY, json[BatchY]);
    fillCommits(m_piLow, json[PiLow]);
    fillCommits(m_piMid, json[PiMid]);
    fillCommits(m_piHigh, json[PiHigh]);

    return true;
}

in::SplinePolynom ProverProof::correctGates(const in::SplittedT_t &t, const CircutParams::SParams_t &SParams) const
{
    const auto &left = t.left.toSplinePolynom();
    const auto &right = t.right.toSplinePolynom();

    auto funcF = in::SplinePolynom{};

    for(const auto &[operation, dots] : SParams.opsFromS) {

        auto isOperation = in::InterpolationPolynom(dots).toSplinePolynom();

        /*Операции через InterpolationPolynom*/
        switch(operation) {
        case BaseGateType::Sum: {
            funcF += (left + right) * isOperation;
            break;
        }
        case BaseGateType::Product: {
            funcF += (left * right) * isOperation;
            break;
        }
        }

    }

    return funcF - t.result.toSplinePolynom();
}

ProverProof::WResult_t ProverProof::correctPermulations(const in::witnesses_t &witnesses, const in::SplinePolynom &num, const in::SplinePolynom &den) const
{
    in::dots_t WDots, WDotsNext;
    WDots.reserve(witnesses.size());
    WDotsNext.reserve(witnesses.size());

    in::Y_t currN, currD, currW = 1;

    WDots.push_back({witnesses.front(), 1});
    for (auto it = witnesses.begin(); it != std::prev(witnesses.end()); it++) {
        currN = num(*it);
        currD = den(*it);

        currW *= currN / currD;
        WDots.push_back({*(it + 1), currW});
        WDotsNext.push_back({*it, currW});
    }
    WDotsNext.push_back({witnesses.back(), 1});

    return {in::InterpolationPolynom(WDots).toSplinePolynom(),
                in::InterpolationPolynom(WDotsNext).toSplinePolynom()};
}

in::SplinePolynom ProverProof:: correctСontinuity(in::SplinePolynom poly) const
{
    std::size_t segmetsCount = poly.segmentsCount();

    for(std::size_t i = 1; i < segmetsCount; i++) {
        poly.at(i - 1).second -= poly.at(i).second;
    }

    poly.at(segmetsCount - 1).second = in::CanonicPolynom::Zero();

    return poly;
}

in::hash_t ProverProof::rWitnessLeaf(int i, std::array<in::value_t, 12> alphaPow) const
{
    const auto &bc = m_batchCommits[i];

    in::commit_t comm = bc.T * alphaPow[0] + bc.WT * alphaPow[1] + bc.WI * alphaPow[2] +
                        bc.W * alphaPow[3] + bc.WNext * alphaPow[4] + bc.G * alphaPow[5] +
                        bc.WitnessZ * alphaPow[6] + bc.C * alphaPow[7] + bc.QP * alphaPow[8] +
                        bc.QG * alphaPow[9] + bc.QC * alphaPow[10] + bc.PartitionZ * alphaPow[11];

    return in::hash(toString(comm));
}

// O(n)
MerkleTree ProverProof::makeTree(in::SplinePolynom &polynom, ProverParams &pp) const
{
    using preLeaf_t = std::vector<std::string>;

    std::size_t segmentsCount = polynom.segmentsCount();

    preLeaf_t preLeafs;
    preLeafs.resize(segmentsCount);

    auto func = [&pp, &polynom, &preLeafs](in::RangeIterator it, in::RangeIterator end)
    {
        for(; it != end; ++it) {
            preLeafs[*it] = toString(polynom.commit(pp, *it));
        }
    };

    in::Threader<in::RangeIterator, preLeaf_t> t({0}, {segmentsCount}, func);

    t();

    return MerkleTree{preLeafs};
}

}
