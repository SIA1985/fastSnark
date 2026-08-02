#ifndef THREADER_H
#define THREADER_H

#include <future>
#include <execution>

#include "params.h"


namespace in {

class RangeIterator
{
public:
    using iterator_concept = std::random_access_iterator_tag;
    using iterator_category = std::random_access_iterator_tag;
    using value_type        = std::size_t;
    using difference_type   = std::ptrdiff_t;
    using pointer           = const std::size_t*;
    using reference         = std::size_t;

    RangeIterator() = default;
    RangeIterator(std::size_t begin, std::size_t step = 1);

    std::size_t step() const;

    RangeIterator &operator++();
    RangeIterator &operator--();
    RangeIterator &operator+=(difference_type n);
    RangeIterator &operator-=(difference_type n);

    bool operator!=(const RangeIterator &other) const;
    bool operator==(const RangeIterator &other) const;
    bool operator<(const RangeIterator &other) const;

    std::size_t operator*() const;

    friend RangeIterator operator+(RangeIterator it, difference_type n) { it += n; return it; }
    friend RangeIterator operator+(difference_type n, RangeIterator it) { it += n; return it; }
    friend RangeIterator operator-(RangeIterator it, difference_type n) { it -= n; return it; }
    friend difference_type operator-(const RangeIterator& a, const RangeIterator& b) { return a.m_current - b.m_current; }

private:
    std::size_t m_current, m_step;
};

#define sortMT(b, e)            snrk::MultiThreading ? std::sort(std::execution::par, b, e) : \
                                                       std::sort(std::execution::seq, b, e)
#define copyMT(b, e, to)        snrk::MultiThreading ? std::copy(std::execution::par, b, e, to) : \
                                                       std::copy(std::execution::seq, b, e, to)
#define transformMT(b, e, to, pred)     snrk::MultiThreading ? std::transform(std::execution::par, b, e, to, pred) : \
                                                               std::transform(std::execution::seq, b, e, to, pred)
#define forEachMT(b, e, pred)   snrk::MultiThreading ? std::for_each(std::execution::par, b, e, pred) : \
                                                       std::for_each(std::execution::seq, b, e, pred)

template <typename iterator_t, typename container_t>
class Threader
{
public:
    using rangeFunc_t = std::function<void(iterator_t, iterator_t, container_t&)>;
    using containers_t = std::vector<container_t>;

    using pair_t = struct{iterator_t begin, end;};
    using pairFunc_t = std::function<void(pair_t, pair_t, container_t&)>;

    using fillFunc_t = std::function<void(iterator_t, iterator_t)>;

    enum Bound : char {
        NoBound,     /*Без захвата границ чанков*/
        WithBound,   /*С захватом границ чанков (для полиномов)*/
    };

    enum Type : char {
        Range,
        TwoIterators,
        Fill,
    };

    Threader(iterator_t begin, iterator_t end, rangeFunc_t func)
        : m_type{Range}
        , m_begin{begin}
        , m_end{end}
        , m_rangeFunc{func}
    {
    }

    Threader(pair_t a, pair_t b, pairFunc_t func)
        : m_type{TwoIterators}
        , m_a{a}
        , m_b{b}
        , m_pairFunc{func}
    {
    }

    Threader(iterator_t begin, iterator_t end, fillFunc_t func)
        : m_type{Fill}
        , m_begin{begin}
        , m_end{end}
        , m_fillFunc{func}
    {
    }

    std::optional<containers_t> operator()(Bound bound = NoBound, std::size_t batchMultiplicity = 1)
    {
        using split_t = std::function<std::pair<std::vector<pair_t>, int>(iterator_t, iterator_t)>;

        std::vector<std::future<void>> futures;
        containers_t containers;

        split_t splitPayload = [batchMultiplicity] (iterator_t begin, iterator_t end) -> std::pair<std::vector<pair_t>, int>
        {
            std::vector<pair_t> splittedPayload;

            std::size_t size = std::distance(begin, end);
            int threadCount = std::thread::hardware_concurrency();

            if (size < static_cast<std::size_t>(threadCount)) {
                threadCount = size > 0 ? static_cast<int>(size) : 1;
            }

            splittedPayload.reserve(threadCount);
            iterator_t it = begin;

            std::size_t baseChunkSize = size / threadCount;
            std::size_t remainder = size % threadCount;

            for (int i = 0; i < threadCount; ++i) {
                auto localBegin = it;

                if (localBegin == end) {
                    splittedPayload.push_back({end, end});
                    continue;
                }

                std::size_t currentChunkSize = baseChunkSize + (i < remainder ? 1 : 0);

                std::advance(it, currentChunkSize);

                splittedPayload.push_back({localBegin, it});
            }

            return {splittedPayload, threadCount};
        };

        split_t splitPayloadWithBound = [batchMultiplicity] (iterator_t begin, iterator_t end) -> std::pair<std::vector<pair_t>, int>
        {
            std::vector<pair_t> splittedPayload;

            auto size = std::distance(begin, end);
            int threadCount = std::thread::hardware_concurrency();

            std::size_t totalSegments = 0;
            if (batchMultiplicity <= 1) {
                totalSegments = size;
            } else {
                totalSegments = (size - 1) / (batchMultiplicity - 1);
            }

            if (totalSegments < static_cast<std::size_t>(threadCount)) {
                threadCount = totalSegments > 0 ? static_cast<int>(totalSegments) : 1;
            }

            splittedPayload.reserve(threadCount);
            iterator_t it = begin;

            std::size_t baseSegmentsPerThread = totalSegments / threadCount;
            std::size_t remainderSegments = totalSegments % threadCount;

            for (int i = 0; i < threadCount; ++i) {
                auto localBegin = it;

                if (i == threadCount - 1) {
                    splittedPayload.push_back({localBegin, end});
                    break;
                }

                std::size_t segmentsForThisThread = baseSegmentsPerThread + (i < remainderSegments ? 1 : 0);

                if (segmentsForThisThread == 0) {
                    threadCount = i;
                    break;
                }

                std::size_t pointsToRead = 0;
                if (batchMultiplicity <= 1) {
                    pointsToRead = segmentsForThisThread;
                } else {
                    pointsToRead = (segmentsForThisThread * (batchMultiplicity - 1)) + 1;
                }

                std::size_t remainingPoints = std::distance(it, end);

                if (remainingPoints <= pointsToRead) {
                    splittedPayload.push_back({localBegin, end});
                    threadCount = i + 1;
                    break;
                }

                std::advance(it, pointsToRead);

                splittedPayload.push_back({localBegin, it});

                std::advance(it, -1);
            }

            return {splittedPayload, threadCount};

        };

        switch (m_type) {
        case Range: {
            if (!snrk::MultiThreading) {
                container_t result;

                m_rangeFunc(m_begin, m_end, result);

                return {containers_t{result}};
            }

            const auto &[splittedPayload, threadCount] =  (bound == NoBound) ? splitPayload(m_begin, m_end) : splitPayloadWithBound(m_begin, m_end);
            containers.resize(threadCount);

            for(int i = 0; i < threadCount; i++) {
                futures.push_back(std::async(m_rangeFunc, splittedPayload[i].begin, splittedPayload[i].end, std::ref(containers[i])));
            }

            break;
        }
        case TwoIterators: {
            if (!snrk::MultiThreading) {
                container_t result;

                m_pairFunc(m_a, m_b, result);

                return {containers_t{result}};
            }

            split_t &split = (bound == NoBound) ? splitPayload : splitPayloadWithBound;

            const auto &[aSplittedPayload, threadCount] = split(m_a.begin, m_a.end);
            containers.resize(threadCount);

            auto bSplittedPayload = split(m_b.begin, m_b.end).first;

            for(int i = 0; i < threadCount; i++) {
                futures.push_back(std::async(m_pairFunc, aSplittedPayload[i], bSplittedPayload[i], std::ref(containers[i])));
            }

            break;
        }

        case Fill: {
            if (!snrk::MultiThreading) {
                container_t result;

                m_fillFunc(m_begin, m_end);

                return std::nullopt;
            }

            const auto &[splittedPayload, threadCount] =  (bound == NoBound) ? splitPayload(m_begin, m_end) : splitPayloadWithBound(m_begin, m_end);

            for(int i = 0; i < threadCount; i++) {
                futures.push_back(std::async(m_fillFunc, splittedPayload[i].begin, splittedPayload[i].end));
            }

            break;
        }

        }

        for(auto &future : futures) {
            future.get();
        }

        switch (m_type) {
        case Fill: {
            return std::nullopt;
        }
        default: {
            return {containers};
        }
        }

    }

private:
    Type m_type;

    iterator_t m_begin, m_end;

    rangeFunc_t m_rangeFunc;
    fillFunc_t m_fillFunc;

    pair_t m_a, m_b;
    pairFunc_t m_pairFunc;

};

}

#endif // THREADER_H
