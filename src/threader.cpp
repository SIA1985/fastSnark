#include "threader.h"
#include "polynom.h"

namespace in {

RangeIterator::RangeIterator(std::size_t begin, std::size_t step)
    : m_current{begin}
    , m_step{step}
{
}

std::size_t RangeIterator::step() const
{
    return m_step;
}

RangeIterator &RangeIterator::operator++()
{
    m_current += m_step;

    return *this;
}

RangeIterator &RangeIterator::operator--()
{
    m_current -= m_step;

    return *this;
}

RangeIterator &RangeIterator::operator+=(difference_type n)
{
    m_current += n * m_step;

    return *this;
}

RangeIterator &RangeIterator::operator-=(difference_type n)
{
    m_current -= n * m_step;

    return *this;
}

bool RangeIterator::operator!=(const RangeIterator &other) const
{
    return m_current != other.m_current;
}

bool RangeIterator::operator==(const RangeIterator &other) const
{
    return m_current == other.m_current;
}

bool RangeIterator::operator<(const RangeIterator &other) const
{
    return m_current < other.m_current;
}

std::size_t RangeIterator::operator*() const
{
    return m_current;
}

}
