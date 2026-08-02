#ifndef POLYNOM_H
#define POLYNOM_H

#include "types.h"

#include "params.h"
#include "threader.h"

namespace in {

using coefs_t = std::vector<value_t>;
using roots_t = xs_t;

class Polynom
{
public:
    Polynom() = default;
    virtual ~Polynom() = default;

    virtual Y_t operator()(X_t x) const = 0;

    virtual commit_t commit(snrk::ProverParams &pp) const;
};

class CanonicPolynom : public Polynom
{
public:
    using devideResult_t = std::pair<CanonicPolynom, CanonicPolynom>;


    CanonicPolynom() = default;

    CanonicPolynom(coefs_t coefs);

    CanonicPolynom(std::size_t n);

    static CanonicPolynom Zero();

    static coefs_t coefsFromRoots(roots_t roots);

    virtual Y_t operator()(X_t x) const override;

    CanonicPolynom operator/(CanonicPolynom &other) const;

    devideResult_t tryDevide(CanonicPolynom &other) const;

    CanonicPolynom operator+(const CanonicPolynom &other) const;

    CanonicPolynom operator+(const value_t value) const;

    CanonicPolynom operator-(const CanonicPolynom &other) const;

    CanonicPolynom operator-(const value_t value) const;

    CanonicPolynom operator*(const CanonicPolynom &other) const;

    CanonicPolynom operator*(const value_t value) const;

    void operator+=(const CanonicPolynom &other);

    void operator+=(const value_t value);

    void operator*=(const CanonicPolynom &other);

    void operator*=(const value_t value);

    void operator-=(const value_t value);

    void operator-=(const CanonicPolynom &other);

    CanonicPolynom operator()(const CanonicPolynom &other) const;

    value_t &operator[](std::size_t i);
    const value_t operator[](std::size_t i) const;

    bool isZero() const;

    virtual commit_t commit(snrk::ProverParams &pp) const override;

protected:
    static CanonicPolynom buildPolynomialRecursive(const roots_t& roots, roots_t::const_iterator start, roots_t::const_iterator end);

    /*x0, x1 .. xn*/
    coefs_t m_coefs;

    friend class SplinePolynom;
};

/* [left, right] */
class Range
{
public:
    enum pos_t : int
    {
        left        = 1 << 0,
        right       = 1 << 1,
        crossed     = 1 << 2,
        equal       = 1 << 3,
        inside      = 1 << 4,
        outside     = 1 << 5,
    };

    //todo:
    Range() = default;

    Range(X_t left, X_t right);

    bool inRangeStrict(X_t x) const;

    bool inRange(X_t x) const;

    pos_t isCrossStrict(const Range &other) const;

    Range crossByStrict(const Range &other) const;

    X_t leftBound() const;
    X_t rightBound() const;

    static Range fromUnsorted(X_t a, X_t b);

private:
    X_t m_left;
    X_t m_right;
};

bool operator<(const Range &a, const Range &b);
bool operator==(const Range &a, const Range &b);
bool operator<=(const Range &a, const Range &b);
std::ostream &operator<<(std::ostream &out, const Range &r);

/*Непрерывный диапазон с одинаковым шагом*/
template<typename T>
class RangeMap
{
public:
    using pair_t = std::pair<Range, T>;
    using data_t = std::vector<pair_t>;

    using const_iterator = typename data_t::const_iterator;
    using iterator = typename data_t::iterator;

    RangeMap() = default;

//    RangeMap(std::size_t)

    RangeMap(const_iterator begin, const_iterator end)
    {
        m_data.resize(std::distance(begin, end));

        copyMT(begin, end, m_data.begin());
    }

    const_iterator find(X_t x) const
    {
        assert(!m_data.empty());

        auto first = m_data.begin();
        if (x < first->first.leftBound()) {
           return first;
        }

        auto last = std::prev(m_data.end());
        if (last->first.rightBound() < x) {
           return last;
        }

        std::size_t segmentIndex = getUL(x - first->first.leftBound()) / (snrk::SplinePartition - 1);

        if (segmentIndex >= m_data.size()) {
           segmentIndex = m_data.size() - 1;
        }

        return first + segmentIndex;
    }

    T operator[](X_t x) const
    {
        return find(x)->second;
    }

    /*Непрерывное заполнение*/
    void insert(const Range &range, const T &polynom)
    {
         if (m_data.size() == 0 ||
             range.leftBound() == std::prev(cend())->first.rightBound()) {
             m_data.push_back(pair_t(range, polynom));
             return;
         }

         if (range.rightBound() == cbegin()->first.leftBound()) {
             m_data.insert(m_data.begin(), pair_t(range, polynom));
             return;
         }

         assert(false);
    }

    void insert(const pair_t &pair)
    {
        insert(pair.first, pair.second);
    }

    void merge(RangeMap &other)
    {
        iterator it;

        if (size() == 0) {
            it = m_data.begin();
        } else
        if (other.size() == 0) {
            return;
        } else
        if (other.begin()->first.leftBound() == std::prev(end())->first.rightBound()) {
            it = m_data.end();
            if (std::size_t expextedSize = other.size() + size(); m_data.capacity() >= expextedSize) {
                m_data.resize(expextedSize);
                copyMT(other.begin(), other.end(), it);
                return;
            }
        } else
        if (std::prev(other.end())->first.rightBound() == begin()->first.leftBound()) {
            it = m_data.begin();
        } else {
            assert(false);
        }

        m_data.insert(it, other.begin(), other.end());

        other.m_data.clear();
    }

    void reserve(std::size_t n)
    {
        m_data.reserve(n);
    }

    const_iterator cbegin() const
    {
        return m_data.cbegin();
    }

    const_iterator cend() const
    {
        return m_data.cend();
    }

    iterator begin()
    {
        return m_data.begin();
    }

    iterator end()
    {
        return m_data.end();
    }

    std::size_t size() const
    {
        return m_data.size();
    }

    const pair_t &atSegment(std::size_t index) const
    {
        return m_data.at(index);
    }

    pair_t &atSegment(std::size_t index)
    {
        return m_data[index];
    }

private:
    data_t m_data;
};

class SplinePolynom : public Polynom
{
public:
    using map_t = RangeMap<CanonicPolynom>;

    using SplineConstIterator = map_t::const_iterator;

    using segment_t = std::size_t;


    SplinePolynom() = default;

    SplinePolynom(const RangeMap<CanonicPolynom> &map);

    /*O(n)*/
    SplinePolynom(dots_t dots, bool fromInterpolation = true);

    /*O(1)*/
    virtual Y_t operator()(X_t x) const override;

    CanonicPolynom operator[](X_t x) const;

    SplinePolynom operator+(const SplinePolynom &other) const;

    SplinePolynom operator+(const value_t value) const;

    SplinePolynom operator-(const SplinePolynom &other) const;

    SplinePolynom operator-(const value_t value) const;

    SplinePolynom operator*(const SplinePolynom &other) const;

    SplinePolynom operator*(const value_t value) const;

    SplinePolynom operator/(SplinePolynom &other) const;

    SplinePolynom operator/(CanonicPolynom &other) const;

    void operator+=(const SplinePolynom &other);

    std::size_t segmentsCount() const;

    segment_t segment(X_t x) const;

    map_t::pair_t &at(segment_t segment);

    // O(n)
    commit_t commit(snrk::ProverParams &pp, segment_t segmentIndex) const;

protected:
    using operatorPred_t = std::function<CanonicPolynom(RangeMap<CanonicPolynom>::const_iterator it,
                                              RangeMap<CanonicPolynom>::const_iterator itOther)>;
    map_t operatorPrivate(const SplinePolynom &other, operatorPred_t pred) const;

    map_t m_map;
};

class InterpolationPolynom : public Polynom
{
public:
    using const_iterator = dots_t::const_iterator;


    InterpolationPolynom() = default;

    InterpolationPolynom(const dots_t &dots);

    /*O(n)*/
    virtual Y_t operator()(X_t x) const override;

    /*O(1)*/
    Y_t operator()(witness_t w) const;

    void operator+=(value_t v);

    void operator*=(value_t v);

    void operator+=(const InterpolationPolynom &other);

    void operator*=(const InterpolationPolynom &other);

    void operator+=(const SplinePolynom &other);

    void operator-=(const InterpolationPolynom &other);

    CanonicPolynom toCanonicPolynom() const;

    SplinePolynom toSplinePolynom() const;

protected:
    dots_t m_dots;
};

class ZeroWitnessPolynom : public CanonicPolynom
{
public:
    ZeroWitnessPolynom(const witnesses_t &xs);

    SplinePolynom toSplinePolynom() const;

    static SplinePolynom makePartitionZeroPolynom(const in::SplinePolynom &WitnessZ);

private:
    witnesses_t m_roots;
};

}
#endif // POLYNOM_H
