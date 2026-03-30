#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace btllib {

struct SerializedMinimizers
{
  std::vector<std::uint8_t> bytes;
  std::size_t count = 0;
};

SerializedMinimizers
minimize_sequence(const std::string& seq,
                  std::size_t k,
                  std::size_t w,
                  std::uint16_t record_idx,
                  std::uint16_t assembly_idx,
                  std::uint8_t is_target);

} // namespace btllib
