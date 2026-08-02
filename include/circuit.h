#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <mutex>

#include "types.h"

namespace snrk {

enum class BaseGateType : char
{
    Sum     = 0,
    Product = 1,
};

/*Так как все приводится к умножению*/
#define IS_SCALED(type) (type == snrk::BaseGateType::Product)

#define FOROPS for(auto operation : {snrk::BaseGateType::Sum, snrk::BaseGateType::Product})

class CircuitValue
{
    using ptr_t = std::shared_ptr<in::value_t>;

public:
    CircuitValue();
    CircuitValue(int data);
    CircuitValue(long long data);
    CircuitValue(double data);
    CircuitValue(unsigned long long data);

    /*Ссылается на один и тот же объект*/
    CircuitValue(const CircuitValue &value);

    CircuitValue& operator=(const CircuitValue&) = delete;

    bool operator<(const CircuitValue &value);

    in::address_t address() const;

    in::value_t value() const;

    static in::value_t scale(in::value_t toScale);

private:
    void fromInt(const in::value_t &integer);

    ptr_t m_data{ptr_t(new in::value_t)};
};
typedef std::vector<CircuitValue> CircutValues;

class BaseGate
{
    using input_t = struct{CircuitValue a; CircuitValue b;};

public:
    BaseGate(BaseGateType type, input_t input, CircuitValue output);

    static const int size;

private:
    BaseGateType m_type{};
    input_t m_input{};
    CircuitValue m_output{};

    friend class CircutParams;
    friend class Circuit;
};

class Circuit
{
    using input_t = std::vector<CircuitValue>;
    using basegates_t = std::vector<BaseGate>;

    using lock_t = std::lock_guard<std::mutex>;

    using gateInput_t = struct{CircuitValue a; CircuitValue b;};

public:
    Circuit(const input_t &input);

    std::size_t size();

    std::size_t degree();


    void addGate(GateType type, const gateInput_t &input, const CircuitValue &output);

private:
    std::mutex m_mutex;

    basegates_t m_gates;

    friend class CircutParams;
};

}

#endif // CIRCUIT_H
