#include "polynom.h"
#include "threader.h"

#include <csignal>

namespace in {

commit_t Polynom::commit(snrk::ProverParams &pp) const
{
    /*Полином не поддерживает обязательство*/
    assert(false);
    return {};
}

CanonicPolynom::CanonicPolynom(coefs_t coefs)
    : m_coefs{coefs}
{
}

CanonicPolynom::CanonicPolynom(std::size_t n)
{
    m_coefs.resize(n, 0.0);
}

CanonicPolynom CanonicPolynom::Zero()
{
    return CanonicPolynom(std::size_t{1});
}

CanonicPolynom CanonicPolynom::ZeroPolynom()
{
    roots_t roots;
    roots.reserve(snrk::SplinePartition);

    for(std::size_t r = snrk::wStart; r <= snrk::SplinePartition; r += snrk::wStep) {
        roots.push_back(r);
    }

    return CanonicPolynom(CanonicPolynom::coefsFromRoots(roots));
}

CanonicPolynom CanonicPolynom::buildPolynomialRecursive(const roots_t& roots, roots_t::const_iterator start, roots_t::const_iterator end)
{
    if (start == end) {
        return {size_t(1)};
    }

    if (std::next(start) == end) {
        CanonicPolynom binomial(2);
        binomial[0] = -(*start);
        binomial[1] = 1.0;
        return binomial;
    }

    auto mid = start;
    std::advance(mid, std::distance(start, end) / 2);

    CanonicPolynom leftPoly = buildPolynomialRecursive(roots, start, mid);

    CanonicPolynom rightPoly = buildPolynomialRecursive(roots, mid, end);

    return leftPoly * rightPoly;
}

coefs_t CanonicPolynom::coefsFromRoots(roots_t roots)
{
    if (roots.empty()) {
        return {value_t(1)};
    }

    CanonicPolynom resultPoly = buildPolynomialRecursive(roots, roots.begin(), roots.end());
    return resultPoly.m_coefs;
}

Y_t CanonicPolynom::operator()(X_t x) const
{
    Y_t y = 0;
    X_t xPow = 1;

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        y += xPow * m_coefs[i];
        xPow *= x;
    }

    return y;
}

bool CanonicPolynom::isZero() const
{
    bool result = true;
    for(const auto &coef : m_coefs) {
        result = result && (coef == value_t(0));

        if(!result) {
            break;
        }
    }

    return result;
}

commit_t CanonicPolynom::commit(snrk::ProverParams &pp) const
{
    /*pp.keys = [1, tg, t^2g...]*/
    assert(pp.keys.size() >= m_coefs.size());

    commit_t com;
    G1::mulVec(com, pp.keys.data(), m_coefs.data(), m_coefs.size());

    return com;
}

CanonicPolynom CanonicPolynom::operator/(CanonicPolynom &other) const
{
    const auto &[res, rem] = tryDevide(other);

    if (!rem.isZero()) {
        std::raise(SIGFPE);
    }

    return res;
}

CanonicPolynom CanonicPolynom::operator+(const CanonicPolynom &other) const
{
    coefs_t diff;
    size_t max_size = std::max(m_coefs.size(), other.m_coefs.size());
    diff.resize(max_size, 0.0);

    for (size_t i = 0; i < max_size; ++i) {
        value_t c1 = (i < m_coefs.size()) ? m_coefs[i] : value_t(0.0);
        value_t c2 = (i < other.m_coefs.size()) ? other.m_coefs[i] : value_t(0.0);
        diff[i] = c1 + c2;
    }

    return CanonicPolynom(diff);
}

CanonicPolynom CanonicPolynom::operator+(const value_t value) const
{
    assert(m_coefs.size() > 0);

    coefs_t result = m_coefs;
    result[0] += value;

    return CanonicPolynom(result);

}

CanonicPolynom CanonicPolynom::operator-(const CanonicPolynom &other) const
{
    coefs_t diff;
    size_t max_size = std::max(m_coefs.size(), other.m_coefs.size());
    diff.resize(max_size, 0.0);

    for (size_t i = 0; i < max_size; ++i) {
        value_t c1 = (i < m_coefs.size()) ? m_coefs[i] : value_t(0.0);
        value_t c2 = (i < other.m_coefs.size()) ? other.m_coefs[i] : value_t(0.0);
        diff[i] = c1 - c2;
    }

    return CanonicPolynom(diff);
}

CanonicPolynom CanonicPolynom::operator-(const value_t value) const
{
    return this->operator+(-value);
}

CanonicPolynom CanonicPolynom::operator*(const CanonicPolynom &other) const
{
    std::size_t newDegree = m_coefs.size() + other.m_coefs.size() - 1;
    CanonicPolynom result(newDegree);

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        for(std::size_t j = 0; j < other.m_coefs.size(); j++) {
            result[i + j] += m_coefs[i] * other.m_coefs[j];
        }
    }
    return result;
}

CanonicPolynom CanonicPolynom::operator*(const value_t value) const
{
    coefs_t result(m_coefs.size(), 0);

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        result[i] = m_coefs[i] * value;
    }

    return CanonicPolynom(result);
}

void CanonicPolynom::operator+=(const CanonicPolynom &other)
{
    std::size_t nOther = other.m_coefs.size();
    m_coefs.resize(std::max(nOther, m_coefs.size()), 0.0);

    for(std::size_t i = 0; i < nOther; i++) {
        m_coefs[i] += other.m_coefs[i];
    }
}

void CanonicPolynom::operator+=(const value_t value)
{
    assert(m_coefs.size() > 0);

    m_coefs[0] += value;
}

void CanonicPolynom::operator*=(const CanonicPolynom &other)
{
    m_coefs = std::move((*this * other).m_coefs);
}

void CanonicPolynom::operator-=(const value_t value)
{
    assert(m_coefs.size() > 0);

    m_coefs[0] -= value;
}

void CanonicPolynom::operator-=(const CanonicPolynom &other)
{
    std::size_t nOther = other.m_coefs.size();
    m_coefs.resize(std::max(nOther, m_coefs.size()), 0.0);

    for(std::size_t i = 0; i < nOther; i++) {
        m_coefs[i] -= other.m_coefs[i];
    }
}

void CanonicPolynom::operator*=(const value_t value)
{
    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        m_coefs[i] *= value;
    }
}

CanonicPolynom CanonicPolynom::operator()(const CanonicPolynom &other) const
{
    if (m_coefs.empty()) {
        return {};
    }

    CanonicPolynom result = CanonicPolynom(coefs_t{m_coefs.back()});

    for (int i = m_coefs.size() - 2; i >= 0; --i) {
        result *= other;
        result += CanonicPolynom(coefs_t{m_coefs[i]});
    }

    return result;
}

value_t &CanonicPolynom::operator[](std::size_t i)
{
    return m_coefs[i];
}

const value_t CanonicPolynom::operator[](std::size_t i) const
{
    return m_coefs.size() < i ? value_t(0.0) : m_coefs.at(i);
}

CanonicPolynom::devideResult_t CanonicPolynom::tryDevide(CanonicPolynom &other) const
{
    assert(!(other.m_coefs.empty() ||
            (other.m_coefs.size() == 1 && other.m_coefs[0] == value_t(0))));

    std::size_t n = m_coefs.size() - 1;
    std::size_t nOther = other.m_coefs.size() - 1;

    if (n < nOther) {
       CanonicPolynom res = CanonicPolynom::Zero();
       CanonicPolynom rem = *this;
       return {res, rem};
    }

    std::size_t countRes = n - nOther + 1;
    CanonicPolynom res(countRes);

    coefs_t remainderCoefs(m_coefs.begin(), m_coefs.end());

    int degreeRes = n - nOther;
    for (int i = degreeRes; i >= 0; --i) {
       value_t factor = remainderCoefs[i + nOther] / other.m_coefs[nOther];
       res[i] = factor;

       for (int j = 0; j <= nOther; ++j) {
           remainderCoefs[i + j] -= factor * other.m_coefs[j];
       }
    }

    std::size_t countRem = (nOther == 0) ? 1 : nOther;
    CanonicPolynom rem(countRem);

    if (nOther > 0) {
       for (std::size_t j = 0; j < countRem; ++j) {
           rem[j] = remainderCoefs[j];
       }
    }

    return {res, rem};
}

InterpolationPolynom::InterpolationPolynom(const dots_t &dots)
    : m_dots{dots}
{
    //todo: убрать лишние соритровки
    sortMT(m_dots.begin(), m_dots.end());
}

Y_t InterpolationPolynom::operator()(X_t x) const
{
//    auto l = [this](std::size_t i, const X_t &x)
//    {
//        X_t li = 1;
//        for(std::size_t j = 0; j < m_dots.size(); j++) {
//            if (i == j) {
//                continue;
//            }

//            li *= ((x - m_dots[j].x) / (m_dots[i].x - m_dots[j].x));
//        }

//        return li;
//    };

//    Y_t y = 0;

//    for(std::size_t i = 0; i < m_dots.size(); i++) {
//        y += (l(i, x) * m_dots[i].y);
//    }

//    return y;
    assert(false);
    return 0;
}

InterpolationPolynom InterpolationPolynom::operator+(const InterpolationPolynom &other) const
{
    return operatorPrivate(other, []OPERATORHEADER
    {
        return a + b;
    });
}

InterpolationPolynom InterpolationPolynom::operator*(const InterpolationPolynom &other) const
{
    return operatorPrivate(other, []OPERATORHEADER
    {
        return a * b;
    });
}

Y_t InterpolationPolynom::operator()(witness_t w) const
{
    std::size_t index = (w - snrk::wStart) / snrk::wStep;

    return m_dots.at(index).y;
}

void InterpolationPolynom::operator+=(value_t v)
{
    for(auto &dot : m_dots) {
        dot.y += v;
    }
}

void InterpolationPolynom::operator*=(const InterpolationPolynom &other)
{
    assert(m_dots.size() == other.m_dots.size());

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        m_dots[i].y *= other.m_dots[i].y;
    }
}

void InterpolationPolynom::operator+=(const SplinePolynom &other)
{
    for(auto &dot : m_dots) {
        dot.y += other(dot.x);
    }
}

void InterpolationPolynom::operator*=(value_t v)
{
    for(auto &dot : m_dots) {
        dot.y *= v;
    }
}

void InterpolationPolynom::operator+=(const InterpolationPolynom &other)
{
    assert(m_dots.size() == other.m_dots.size());

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        m_dots[i].y += other.m_dots[i].y;
    }
}

void InterpolationPolynom::operator-=(const InterpolationPolynom &other)
{
    assert(m_dots.size() == other.m_dots.size());

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        m_dots[i].y -= other.m_dots[i].y;
    }
}

CanonicPolynom InterpolationPolynom::toCanonicPolynom() const
{
    if (m_dots.empty()) {
        return {};
    }

    int n = m_dots.size();

    //Перевод в локальные координаты
    dots_t dots(m_dots.size());
    for(std::size_t i = 0, x = snrk::wStart; i < m_dots.size(); i++, x += snrk::wStep) {
        dots[i] = {x, m_dots[i].y};
    }

    // Шаг 1: Вычисление разделённых разностей
    std::vector<values_t> divDiff(n, values_t(n));
    for (int i = 0; i < n; ++i) {
        divDiff[i][0] = dots[i].y;
    }

    for (int j = 1; j < n; ++j) {
        for (int i = j; i < n; ++i) {
            divDiff[i][j] = (divDiff[i][j-1] - divDiff[i-1][j-1]) / (dots[i].x - dots[i-j].x);
        }
    }

    // Шаг 2: Получение коэффициентов полинома Ньютона (диагональные элементы)
    coefs_t newtonCoeffs(n);
    for (int i = 0; i < n; ++i) {
        newtonCoeffs[i] = divDiff[i][i];
    }

    // Шаг 3: Преобразование в каноническую форму
    coefs_t canonicalCoeffs(n, 0.0);

    canonicalCoeffs[0] = newtonCoeffs[0];

    coefs_t newtonBasisPoly(n, 0.0);
    newtonBasisPoly[0] = 1.0;

    for (int i = 1; i < n; ++i) {
        coefs_t nextBasisPoly(n, 0.0);

        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k+1] += newtonBasisPoly[k];
        }

        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k] -= newtonBasisPoly[k] * dots[i-1].x;
        }
        newtonBasisPoly = nextBasisPoly;

        for (int j = 0; j <= i; ++j) {
            canonicalCoeffs[j] += newtonCoeffs[i] * newtonBasisPoly[j];
        }
    }

    return CanonicPolynom(canonicalCoeffs);
}

SplinePolynom InterpolationPolynom::toSplinePolynom() const
{
    return SplinePolynom(m_dots);
}

InterpolationPolynom InterpolationPolynom::operatorPrivate(const InterpolationPolynom &other, std::function<Y_t(const Y_t &, const Y_t &)> pred) const
{
    assert(m_dots.front().x == other.m_dots.front().x);

    std::size_t min = std::min(m_dots.size(), other.m_dots.size()),
                max = std::max(m_dots.size(), other.m_dots.size());

    const dots_t &longest = (m_dots.size() > other.m_dots.size()) ? m_dots : other.m_dots;

    dots_t dots(max);

    for(std::size_t i = 0; i < min; i++) {
        dots[i] = {m_dots[i].x, pred(m_dots[i].y, other.m_dots[i].y)};
    }

    for(std::size_t i = min; i < max; i++) {
        dots[i] = longest[i];
    }

    return InterpolationPolynom(dots);
}

Range::Range(X_t left, X_t right)
    : m_left{left}
    , m_right{right}
{
    assert(right >= left); //!=-1
}

X_t Range::leftBound() const
{
    return m_left;
}

X_t Range::rightBound() const
{
    return m_right;
}

Range Range::fromUnsorted(X_t a, X_t b)
{
    if (a < b) {
        return {a, b};
    }

    return {b, a};
}

bool operator<(const Range &a, const Range &b)
{
    return a.rightBound() < b.leftBound();
}

bool operator==(const Range &a, const Range &b)
{
    return a.leftBound() == b.leftBound() && a.rightBound() == b.rightBound();
}

bool operator<=(const Range &a, const Range &b)
{
    return a < b || a == b;
}

std::ostream &operator<<(std::ostream &out, const Range &r)
{
    return out << "{" << r.leftBound() << ", " << r.rightBound() << "}";
}

SplinePolynom::SplinePolynom(const map_t &map)
    : m_map{map}
{
}

SplinePolynom::SplinePolynom(dots_t dots)
{
    sortMT(dots.begin(), dots.end());

    using iterator = dots_t::const_iterator;
    using threader = Threader<iterator, map_t>;

    auto fromFunc = [](map_t& m, const dots_t &dots)
    {
        m.insert({dots.front().x, dots.back().x}, InterpolationPolynom(dots).toCanonicPolynom());
    };

    auto func = [fromFunc](iterator begin, iterator end, map_t& m) {
        for(auto it = begin; it != end; /*it = (it == end) ? it : std::prev(it)*/) {
            auto start = it;
            if (std::distance(it, end) >= snrk::SplinePartition) {
                std::advance(it, snrk::SplinePartition);
            } else {
                it = end;
            }

            dots_t dots(start, it);

            fromFunc(m, dots);
        }
    };

    threader t(dots.begin(), dots.end(), func);

    auto maps = t(snrk::SplinePartition).value();

    std::size_t size = 0;
    for(const auto &map : maps) {
        size += map.size();
    }

    m_map.reserve(size);
    for(auto &map : maps) {
        m_map.merge(map);
    }
}

Y_t SplinePolynom::operator()(X_t x) const
{
    if (x > std::prev(m_map.cend())->first.rightBound()) {
        return m_map[x](x);
    }

    value_t localX = getUL(x) % (snrk::SplinePartition + 1);
    return m_map[x](localX);
}

CanonicPolynom SplinePolynom::operator[](X_t x) const
{
    return m_map[x];
}

SplinePolynom SplinePolynom::operator/(SplinePolynom &other) const
{
    return SplinePolynom(operatorPrivate(other, []
    (map_t::const_iterator it, map_t::const_iterator itOther)
    {
        CanonicPolynom f1 = it->second;
        CanonicPolynom f2 = itOther->second;

        auto res = f1 / f2;

        return res;
    }));
}

SplinePolynom SplinePolynom::operator/(CanonicPolynom &other) const
{
    using iterator_t = map_t::iterator;
    auto map = m_map;

    Threader<iterator_t, map_t> t(map.begin(), map.end(),
    [&other, b = map.cbegin()](iterator_t it, iterator_t end)
    {
        for(; it != end; it++) {
            it->second = it->second / other;
        }

    });

    t();

    return SplinePolynom(map);
}

SplinePolynom SplinePolynom::operator+(const SplinePolynom &other) const
{
    return SplinePolynom(operatorPrivate(other, []
    (map_t::const_iterator it, map_t::const_iterator itOther)
    {
        return it->second + itOther->second;
    }));
}

SplinePolynom SplinePolynom::operator+(const value_t value) const
{
    using iterator_t = map_t::iterator;
    auto map = m_map;

    Threader<iterator_t, map_t> t(map.begin(), map.end(),
    [value, b = map.cbegin()](iterator_t it, iterator_t end)
    {
        for(; it != end; it++) {
            it->second = it->second + value;
        }

    });

    t();

    return SplinePolynom(map);
}

SplinePolynom SplinePolynom::operator*(const SplinePolynom &other) const
{
    return SplinePolynom(operatorPrivate(other,
    [](map_t::const_iterator it, map_t::const_iterator itOther)
    {
        return it->second * itOther->second;
    }));
}

SplinePolynom SplinePolynom::operator*(const value_t value) const
{
    using iterator_t = map_t::iterator;
    auto map = m_map;

    Threader<iterator_t, map_t> t(map.begin(), map.end(),
    [value, b = map.cbegin()](iterator_t it, iterator_t end)
    {
        for(; it != end; it++) {
            it->second = it->second * value;
        }

    });

    t();

    return SplinePolynom(map);
}

SplinePolynom SplinePolynom::operator-(const SplinePolynom &other) const
{
    return SplinePolynom(operatorPrivate(other,
    [](map_t::const_iterator it, map_t::const_iterator itOther)
    {
        return it->second - itOther->second;
    }));
}

SplinePolynom SplinePolynom::operator-(const value_t value) const
{
    return this->operator+(-value);
}

void SplinePolynom::operator+=(const SplinePolynom &other)
{
    if (m_map.size() == 0) {
        m_map = other.m_map;
        return;
    }

    *this = std::move(*this + other);
}

SplinePolynom::map_t SplinePolynom::operatorPrivate(const SplinePolynom &other, operatorPred_t pred) const
{
    assert(m_map.size() == other.m_map.size());


    using iterator = map_t::const_iterator;
    using threader = Threader<iterator, map_t>;
    using pair_t = threader::pair_t;

    auto func = [&pred](pair_t a, pair_t b, map_t &result) {
        Range currentRange = {0, 0};
        auto it = a.begin;
        auto itOther = b.begin;

        auto end = a.end;
        auto otherEnd = b.end;

        auto itRange = [&it](){return it->first;};
        auto otherRange = [&itOther](){ return itOther->first;};

        auto predCall = [&it, &itOther, &currentRange, &pred, &result]()
        {
            result.insert(currentRange, pred(it, itOther));
        };

        while(it != end || itOther != otherEnd) {
            currentRange = itRange();
            predCall();

            it++;
            itOther++;

        }
    };

    threader t({m_map.cbegin(), m_map.cend()}, {other.m_map.cbegin(), other.m_map.cend()}, func);

    auto maps = t().value();

    std::size_t size = 0;
    for(const auto &map : maps) {
        size += map.size();
    }

    map_t result;
    result.reserve(size);
    for(auto &map : maps) {
        result.merge(map);
    }

    return result;
}

std::size_t SplinePolynom::segmentsCount() const
{
    return m_map.size();
}

SplinePolynom::segment_t SplinePolynom::segment(X_t x) const
{
    return std::distance(m_map.cbegin(), m_map.find(x));
}

SplinePolynom::map_t::pair_t &SplinePolynom::at(segment_t segment)
{
    return m_map.atSegment(segment);
}

commit_t SplinePolynom::commit(snrk::ProverParams &pp, segment_t segmentIndex) const
{
    return m_map.atSegment(segmentIndex).second.commit(pp);
}

//todo: многопоток
CanonicPolynom SplinePolynom::toCanonicPolynom() const
{
    in::coefs_t coefs;
    coefs.reserve(m_map.size() > 0 ? m_map.size() * m_map.cbegin()->second.m_coefs.size() : 0);

    for(std::size_t i = 0; i < m_map.size(); i++) {
        for(auto c : m_map.atSegment(i).second.m_coefs) {
            coefs.push_back(c);
        }
    }

    return CanonicPolynom(coefs);
}

}
