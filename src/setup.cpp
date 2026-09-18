#include "setup.h"
#include "threader.h"


namespace in {
witnesses_t genWitnesses(witness_t start, std::size_t count, in::witness_t wStep)
{
    witnesses_t witnesses;
    for(witness_t i = 0; i < count; i++) {
        witnesses.push_back(start);
        start += wStep;
    }

    return witnesses;
}

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

}

namespace snrk {

GlobalParams::GlobalParams(in::value_t t)
{
    mcl::mapToG1(m_vp.g1, 1);

    in::G1 current_power = m_vp.g1;
    m_pp.keys.push_back(current_power);

    //todo: чтение глобальных параметров
    for (int i = 1; i < 1'000'000; i++) {
        in::G1::mul(current_power, current_power, t);
        m_pp.keys.push_back(current_power);
    }
    /**/

    mcl::mapToG2(m_vp.g2, 1);
    in::G2::mul(m_vp.tG2, m_vp.g2, t);
}

ProverParams GlobalParams::PP() const
{
    return m_pp;
}

VerifierParams GlobalParams::VP() const
{
    return m_vp;
}

//todo:
void to_json(json_t& j, const GlobalParams& gp)
{
//    j["t"] = tG.t;
//    j["G"] = tG.G;
}

//todo:
void from_json(const snrk::json_t& j, GlobalParams& gp)
{
//    tG.t = j.at("t").get<value_t>();
//    tG.G = j.at("G");
}

CircutParams::CircutParams(Circuit &circut, ProverParams &pp)
    : m_pp{pp}
    , m_circutSize{circut.size()}
{
    m_witnesses = in::genWitnesses(wStart, circut.degree(), wStep);

    if (MultiThreading) {
        std::thread tT(&CircutParams::generateT, this, std::ref(circut));
        std::thread tS(&CircutParams::generateS, this, std::ref(circut));
        std::thread tW(&CircutParams::generateW, this, std::ref(circut));

        tT.join();
        tS.join();
        tW.join();
    } else {
        generateT(circut);
        generateS(circut);
        generateW(circut);
    }
}

in::witnesses_t CircutParams::witnesses() const
{
    return m_witnesses;
}

std::size_t CircutParams::witnessesCount() const
{
    return m_witnesses.size();
}

CircutParams::params_t CircutParams::params() const
{
    return {.TParams = {m_T, m_splittedT}, .SParams = {m_opsFromS}, .WParams = {m_WT, m_WI, map}};
}

std::size_t CircutParams::circutSize() const
{
    return m_circutSize;
}

void CircutParams::generateT(Circuit &circut)
{
    in::dots_t dots, leftDots, rightDots, resultDots;

    auto cw = m_witnesses.cbegin();

    for(const auto& gate : circut.m_gates) {
        auto left = gate.m_input.a.value();
        auto right = gate.m_input.b.value();
        auto result = IS_SCALED(gate.m_type) ? CircuitValue::scale(gate.m_output.value()) : gate.m_output.value();

        dots.push_back({*cw, left});
        leftDots.push_back({*cw, left});
        rightDots.push_back({*cw, right});
        resultDots.push_back({*cw, result});
        cw++;

        dots.push_back({*cw, right});
        leftDots.push_back({*cw, left});
        rightDots.push_back({*cw, right});
        resultDots.push_back({*cw, result});
        cw++;

        dots.push_back({*cw, result});
        leftDots.push_back({*cw, left});
        rightDots.push_back({*cw, right});
        resultDots.push_back({*cw, result});
        cw++;
    }

    m_T = in::T_t(dots);
    m_splittedT = {
                    .left = in::InterpolationPolynom(leftDots),
                    .right = in::InterpolationPolynom(rightDots),
                    .result = in::InterpolationPolynom(resultDots)
                  };
}

void CircutParams::generateS(Circuit &circut)
{
    auto cw = m_witnesses.cbegin();
    for(std::size_t i = 0; i < circut.size(); i++) {
        auto currentOperation = circut.m_gates[i].m_type;

        for(int i = 0; i < BaseGate::size; i++) {
            FOROPS {
                m_opsFromS[operation].push_back(
                            operation == currentOperation ? in::dot_t{*cw, 1} : in::dot_t{*cw, 0}
                );
            }

            cw++;
        }
    }
}

void CircutParams::generateW(Circuit &circut)
{
    /*addr -> свидетели*/
    std::unordered_map<std::size_t, std::shared_ptr<cond_t>> duplicates;
    auto insert = [&duplicates](const CircuitValue &key, in::witness_t value)
    {
        auto address = key.address();
        if (duplicates.count(address) == 0) {
            duplicates[address] = std::make_shared<cond_t>();
        }

        duplicates[address]->insert(value);
    };

    auto cw = m_witnesses.cbegin();

    for(const auto &gate : circut.m_gates) {
        insert(gate.m_input.a, *cw++);
        insert(gate.m_input.b, *cw++);
        insert(gate.m_output, *cw++);
    }


    auto circleIterator = [](cond_t::iterator begin, cond_t::iterator it, cond_t::iterator end)
    {
        if (it == end) {
            return begin;
        }

        return it;
    };

    in::dots_t dotsWI;
    in::dots_t dotsWT;
    dotsWI.reserve(duplicates.size());
    dotsWT.reserve(duplicates.size());

    for(const auto &[_, condition] : duplicates) {
        auto begin = condition->begin();
        auto end = condition->end();

        for(auto it = begin; it != end;) {
            auto current_node = *it;
            auto next_node = *circleIterator(begin, ++it, end);

            dotsWI.push_back({current_node, current_node});
            dotsWT.push_back({current_node, next_node});

            map[current_node] = next_node;
        }
    }

    m_WI = in::W_t(dotsWI); //witness -> witness
    m_WT = in::W_t(dotsWT); //witness -> next_this_addres_value_wintess
}

}
