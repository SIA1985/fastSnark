#include "circuit.h"
#include "params.h"


namespace snrk {

CircuitValue::CircuitValue()
{
}

CircuitValue::CircuitValue(int data)
{
    in::value_t integer;
    integer.setStr(std::to_string(data), 10);

    fromInt(integer);
}

CircuitValue::CircuitValue(long long data)
{
    in::value_t integer;
    integer.setStr(std::to_string(data), 10);

    fromInt(integer);
}

CircuitValue::CircuitValue(double data)
{
    double scaleFactor = 1e14;

    long double scaled = static_cast<long double>(data) * scaleFactor;

    scaled = std::round(scaled);

    bool isNeg = scaled < 0;
    if (isNeg) {
        scaled = -scaled;
    }

    std::stringstream ss;
    ss << std::fixed << std::setprecision(0) << scaled;

    m_data->setStr(ss.str(), 10);

    if (isNeg) {
        in::negative(*m_data);
    }
}

CircuitValue::CircuitValue(unsigned long long data)
{
    in::value_t integer;
    integer.setStr(std::to_string(data), 10);

    fromInt(integer);
}

snrk::CircuitValue::CircuitValue(const CircuitValue &value)
{
    m_data = value.m_data;
}

bool CircuitValue::operator<(const CircuitValue &value)
{
    return *m_data < *value.m_data;
}

in::address_t CircuitValue::address() const
{
    return in::address_t(m_data.get());
}

in::value_t CircuitValue::value() const
{
    return *m_data;
}

in::value_t CircuitValue::scale(in::value_t toScale)
{
    in::value_t k;
    k.setStr("100000000000000", 10);

    return toScale * k;
}

void CircuitValue::fromInt(const in::value_t &integer)
{
    *m_data = scale(integer);
}

const int BaseGate::size = 3;
BaseGate::BaseGate(BaseGateType type, input_t input, CircuitValue output)
    : m_type{type}
    , m_input{input}
    , m_output{output}
{

}

Circuit::Circuit(const input_t &input)
{
    for(const auto &i : input) {
        m_gates.push_back(BaseGate{BaseGateType::Sum, {i, 0}, i});
    }
}

std::size_t Circuit::size()
{
    lock_t lg(m_mutex);

    return m_gates.size();
}

std::size_t Circuit::degree()
{
    lock_t lg(m_mutex);

    return 3 * m_gates.size(); //При кастомных гейтах будет меняться "3"
}

void Circuit::addGate(GateType type, const gateInput_t &input, const CircuitValue &output)
{
    if (Gates.count(type) == 0) {
        return;
    }

    lock_t lg(m_mutex);

    for(const auto &gate : Gates[type](input.a, input.b, output)) {
        m_gates.push_back(gate);
    }
}

}
