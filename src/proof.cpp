#include "proof.h"
#include "threader.h"

#include <chrono>
#define SINCE(s) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - s).count()

namespace snrk {

ProverProof::ProverProof(const CircutParams &cp, ProverParams &pp, const std::string &trInit)
{
    const auto &[TParams, SParams, WParams] = cp.params();
    const in::witnesses_t &witnesses = cp.witnesses();

    //Инициализация Транскрипта
    in::Transcript tr(trInit);

    //todo: что-то добавить в транскрипт до? Да нужно как-то защитить WT, WI от подмены
    //из компиляции кода программы и получения обязательств WT, WI, S
    //поэтому до этого добавляем обязательства WT, WI и S
    auto alphaPow = tr.getPowerArray<2>("alpha");
    auto betaPow = tr.getPowerArray<4>("beta");

    in::CanonicPolynom Z = in::CanonicPolynom::ZeroPolynom();

    //Проверка перестановок
    in::SplinePolynom cP = correctPermulations(witnesses, TParams, WParams, alphaPow);

    in::CanonicPolynom P = cP.toCanonicPolynom();
    in::CanonicPolynom QP = P / Z;

    //Проверка вычислений
    in::SplinePolynom cG = correctGates(TParams, SParams);

    in::CanonicPolynom G = cG.toCanonicPolynom();
    in::CanonicPolynom QG = G / Z;

    //Получение точки проверки
    in::value_t r = tr("r");

    //Создание доказательства
    in::CanonicPolynom F = G * betaPow[0] + QG * betaPow[1] + P * betaPow[2] + QP * betaPow[3];

    m_F = F.commit(pp);
    F -= F(r);

    in::CanonicPolynom Rz = in::CanonicPolynom({-r, 1});
    in::CanonicPolynom Q = F / Rz;

    m_Q = Q.commit(pp);

    //Значения проверяемых полиномов
    m_Pr = P(r);
    m_QPr = QP(r);
    m_QGr = QG(r);
    m_Gr = G(r);
}

bool ProverProof::check(const VerifierParams &vp, const std::string &trInit)
{
    in::Transcript tr(trInit);

    tr.getPowerArray<2>("alpha");
    auto betaPow = tr.getPowerArray<4>("beta");

    in::value_t r = tr("r");

    in::value_t Fr = m_Gr * betaPow[0] + m_QGr * betaPow[1] + m_Pr * betaPow[2] + m_QPr * betaPow[3];

    in::GT e1, e2;
    mcl::pairing(e1, m_Q, vp.tG2 - (vp.g2 * r));
    mcl::pairing(e2, m_F - vp.g1 * Fr, vp.g2);
    if (e1 != e2) {
        return false;
    }

    in::CanonicPolynom Z = in::CanonicPolynom::ZeroPolynom();

    if (m_QGr != m_Gr / Z(r)) {
        return false;
    }

    if (m_QPr != m_Pr / Z(r)) {
        return false;
    }

    return true;
}

const std::string Q = "Q", F = "F", Pr = "Pr",
                  QPr = "QPr", Gr = "Gr", QGr = "QGr";

json_t ProverProof::toJson() const
{
    json_t json;

    json[Q] = toString(m_Q);
    json[F] = toString(m_F);

    json[Pr] = toString(m_Pr);
    json[QPr] = toString(m_QPr);
    json[Gr] = toString(m_Gr);
    json[QGr] = toString(m_QGr);

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    if ((!json.contains(Q) || !json[Q].is_string()) &&
        (!json.contains(F) || !json[F].is_string()) &&
        (!json.contains(Pr) || !json[Pr].is_string()) &&
        (!json.contains(QPr) || !json[QPr].is_string()) &&
        (!json.contains(Gr) || !json[Gr].is_string()) &&
        (!json.contains(QGr) || !json[QGr].is_string())) {
        return false;
    }

    m_Q = fromString<in::commit_t>(json[Q].get<std::string>()).value();
    m_F = fromString<in::commit_t>(json[F].get<std::string>()).value();

    m_Pr = fromString<in::value_t>(json[Pr].get<std::string>()).value();
    m_QPr = fromString<in::value_t>(json[QPr].get<std::string>()).value();
    m_Gr = fromString<in::value_t>(json[Gr].get<std::string>()).value();
    m_QGr = fromString<in::value_t>(json[QGr].get<std::string>()).value();

    return true;
}

in::SplinePolynom ProverProof::correctGates(const CircutParams::TParams_t &TParams, const CircutParams::SParams_t &SParams) const
{
    const auto &left = TParams.splittedT.left;
    const auto &right = TParams.splittedT.right;

    auto funcF = in::SplinePolynom{};

    for(const auto &[operation, dots] : SParams.opsFromS) {

        auto isOperation = in::InterpolationPolynom(dots).toSplinePolynom();

        switch(operation) {
        case BaseGateType::Sum: {
            funcF += (left + right).toSplinePolynom() * isOperation;
            break;
        }
        case BaseGateType::Product: {
            funcF += (left * right).toSplinePolynom() * isOperation;
            break;
        }
        }

    }

    return funcF - TParams.splittedT.result.toSplinePolynom();
}

in::SplinePolynom ProverProof::correctPermulations(const in::witnesses_t &witnesses, const CircutParams::TParams_t &TParams,
                                                   const CircutParams::WParams_t &WParams, std::array<in::value_t, 2> trPow) const
{
    std::size_t wSize = witnesses.size();

    in::dots_t numDots, denDots, WDots, WNextDots;
    for(auto ptr : {&numDots, &denDots, &WDots, &WNextDots}) {
        ptr->reserve(wSize);
    }

    in::Y_t currN, currD, currW = 1;
    for (auto w : witnesses) {
        currN = TParams.t(w) + trPow[0] * WParams.wt(w) + trPow[1];
        currD = TParams.t(w) + trPow[0] * WParams.wi(w) + trPow[1];

        in::Y_t prevW = currW;
        currW *= currN / currD;

        numDots.push_back({w, currN});
        denDots.push_back({w, currD});
        WDots.push_back({w, prevW});
        WNextDots.push_back({w, currW});
    }

    const auto &num = in::InterpolationPolynom(numDots).toSplinePolynom();
    const auto &den = in::InterpolationPolynom(denDots).toSplinePolynom();
    const auto &W = in::InterpolationPolynom(WDots).toSplinePolynom();
    const auto &WNext = in::InterpolationPolynom(WNextDots).toSplinePolynom();

    return WNext * den - W * num;
}

// O(n)
MerkleTree ProverProof::makeTree(const in::SplinePolynom &polynom, ProverParams &pp, in::Transcript &tr) const
{
    using preLeaf_t = in::commits_t;

    std::size_t segmentsCount = polynom.segmentsCount();

    preLeaf_t preLeafs;
    preLeafs.resize(segmentsCount);

    auto func = [&pp, &polynom, &preLeafs](in::RangeIterator it, in::RangeIterator end)
    {
        for(; it != end; ++it) {
            preLeafs[*it] = polynom.commit(pp, *it);
        }
    };

    in::Threader<in::RangeIterator, preLeaf_t> t({0}, {segmentsCount}, func);

    t();

    return MerkleTree{preLeafs, tr};
}

}
