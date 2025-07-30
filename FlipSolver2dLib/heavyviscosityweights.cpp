#include "heavyviscosityweights.h"

void HeavyViscosityWeights::multiply(const std::vector<double> &in, std::vector<double> &out) const
{
    ASSERT(in.size() == out.size());

    if (m_data.empty()) {
        out = in;
        return;
    }

    std::vector<Range> ranges = ThreadPool::i()->splitRange(in.size());

    // multiplyThread(Range(0,in.size()),Range(0,m_data.size()),in,out);

    // return;

    for (size_t i = 0; i < ranges.size(); i++) {
        ThreadPool::i()->enqueue(&MatrixWeights::multiplyThread,
                                 this,
                                 ranges.at(i),
                                 ranges.at(i),
                                 std::cref(in),
                                 std::ref(out));
    }
    ThreadPool::i()->wait();
}
