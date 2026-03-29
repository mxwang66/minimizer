#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

namespace btllib {

struct Minimizer
{
  uint64_t min_hash = 0;
  uint64_t out_hash = 0;
  std::size_t pos = 0;
  bool forward = false;
};

void
minimize_sequence(const std::string& seq,
                  std::size_t k,
                  std::size_t w,
                  const std::function<void(const Minimizer&)>& emit);

} // namespace btllib
