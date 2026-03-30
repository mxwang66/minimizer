#include "minimizer.hpp"

#include "nthash.hpp"

#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace btllib {

namespace {

struct Minimizer
{
  std::uint64_t min_hash = 0;
  std::uint64_t out_hash = 0;
  std::size_t pos = 0;
};

inline void
append_record(std::vector<std::uint8_t>& out,
              std::uint64_t out_hash,
              std::uint32_t pos,
              std::uint16_t record_idx,
              std::uint16_t assembly_idx,
              std::uint8_t is_target)
{
  const auto old_size = out.size();
  out.resize(old_size + serialized_record_size);
  auto* record = out.data() + old_size;

  std::memcpy(record + 0, &out_hash, sizeof(out_hash));
  std::memcpy(record + 8, &pos, sizeof(pos));
  std::memcpy(record + 12, &record_idx, sizeof(record_idx));
  std::memcpy(record + 14, &assembly_idx, sizeof(assembly_idx));
  std::memcpy(record + 16, &is_target, sizeof(is_target));
}

inline void
calc_minimizer(const std::vector<Minimizer>& hashed_kmers_buffer,
               const Minimizer*& min_current,
               const std::size_t idx,
               ssize_t& min_idx_left,
               ssize_t& min_idx_right,
               ssize_t& min_pos_prev,
               const std::size_t w,
               std::vector<std::uint8_t>& serialized_minimizers,
               std::size_t& minimizer_count,
               const std::uint16_t record_idx,
               const std::uint16_t assembly_idx,
               const std::uint8_t is_target)
{
  min_idx_left = ssize_t(idx + 1 - w);
  min_idx_right = ssize_t(idx + 1);

  const auto& min_left =
    hashed_kmers_buffer[std::size_t(min_idx_left) % hashed_kmers_buffer.size()];
  const auto& min_right = hashed_kmers_buffer[(std::size_t(min_idx_right) - 1) %
                                               hashed_kmers_buffer.size()];

  if (min_current == nullptr || min_current->pos < min_left.pos) {
    min_current = &min_left;
    for (ssize_t i = min_idx_left; i < min_idx_right; i++) {
      const auto& min_i = hashed_kmers_buffer[std::size_t(i) % hashed_kmers_buffer.size()];
      if (min_i.min_hash <= min_current->min_hash) {
        min_current = &min_i;
      }
    }
  } else if (min_right.min_hash <= min_current->min_hash) {
    min_current = &min_right;
  }

  if (ssize_t(min_current->pos) > min_pos_prev &&
      min_current->min_hash != std::numeric_limits<uint64_t>::max()) {
    if (min_current->pos > std::numeric_limits<std::uint32_t>::max()) {
      throw std::runtime_error("minimizer position exceeds uint32 range");
    }
    min_pos_prev = ssize_t(min_current->pos);
    append_record(serialized_minimizers,
                  min_current->out_hash,
                  static_cast<std::uint32_t>(min_current->pos),
                  record_idx,
                  assembly_idx,
                  is_target);
    ++minimizer_count;
  }
}

} // namespace

SerializedMinimizers
minimize_sequence(const std::string& seq,
                  std::size_t k,
                  std::size_t w,
                  std::uint16_t record_idx,
                  std::uint16_t assembly_idx,
                  std::uint8_t is_target)
{
  if ((k > seq.size()) || (w > seq.size() - k + 1)) {
    return {};
  }

  SerializedMinimizers output;
  output.bytes.reserve((2 * (seq.size() - k + 1) / w) * serialized_record_size);

  std::vector<Minimizer> hashed_kmers_buffer(w + 1);
  ssize_t min_idx_left = -1;
  ssize_t min_idx_right = -1;
  ssize_t min_pos_prev = -1;
  const Minimizer* min_current = nullptr;

  std::size_t idx = 0;
  for (btllib::NtHash nh(seq, 2, k); nh.roll(); ++idx) {
    auto& hk = hashed_kmers_buffer[idx % hashed_kmers_buffer.size()];
    hk = Minimizer{ nh.hashes()[0], nh.hashes()[1], nh.get_pos() };

    if (idx + 1 >= w) {
      calc_minimizer(hashed_kmers_buffer,
                     min_current,
                     idx,
                     min_idx_left,
                     min_idx_right,
                     min_pos_prev,
                     w,
                     output.bytes,
                     output.count,
                     record_idx,
                     assembly_idx,
                     is_target);
    }
  }

  return output;
}

} // namespace btllib
