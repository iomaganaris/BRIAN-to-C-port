/**
 * @file synapses.cpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief Source file containing utilities functions
 */

#include <algorithm>
#include <vector>

template <typename A, typename B>
void zip(
        const std::vector<A> &a,
        const std::vector<B> &b,
        std::vector<std::pair<A,B> > &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

/* Write the first and second element of the pairs in the
 * given zipped vector into a and b. (This assumes that
 * the vectors have equal length)
 */
template <typename A, typename B>
void unzip(
        const std::vector<std::pair<A, B> > &zipped,
        std::vector<A> &a,
        std::vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

template <typename A, typename B>
void sort_second(std::vector<A>& a, std::vector<B>& b) {
    /// Zip the vectors together
    std::vector<std::pair<std::string, int> > zipped;
    zip(a, b, zipped);

    /// Sort the vector of pairs
    std::sort(std::begin(zipped), std::end(zipped),
              [&](const auto &a, const auto &b) {
                  return a.second > b.second;
              });

    /// Write the sorted pairs back to the original vectors
    unzip(zipped, a, b);
}