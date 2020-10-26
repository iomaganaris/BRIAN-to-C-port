//
// Created by Magkanaris Ioannis on 16.10.20.
//

#pragma once

#include <vector>

/**
 * Fill the zipped vector with pairs consisting of the corresponding
 * elements of a and b. (This assumes that the vectors have equal length)
 */
template <typename A, typename B>
void zip(const std::vector<A> &a, const std::vector<B> &b, std::vector<std::pair<A,B> > &zipped);

/**
 * Write the first and second element of the pairs in the given zipped
 * vector into a and b. (This assumes that the vectors have equal length)
 */
template <typename A, typename B>
void unzip(const std::vector<std::pair<A, B> > &zipped, std::vector<A> &a,  std::vector<B> &b);

/**
 * Sort vectors \a and \b based on \b
 */
template <typename A, typename B>
void sort_second(std::vector<A>& a, std::vector<B>& b);
