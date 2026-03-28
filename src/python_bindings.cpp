#include "fasta_reader.hpp"
#include "minimizer.hpp"

#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

template<typename T>
void
append_scalar(std::string& out, T value)
{
  static_assert(std::is_trivially_copyable_v<T>);
  char bytes[sizeof(T)];
  std::memcpy(bytes, &value, sizeof(T));
  out.append(bytes, sizeof(T));
}

std::tuple<std::string, py::list, std::vector<std::size_t>, std::vector<std::size_t>>
indexlr_impl(const std::vector<std::string>& assembly_paths,
             std::size_t kmerlen,
             std::size_t windowsize,
             const std::vector<std::size_t>& assembly_indices,
             const std::vector<bool>& is_targets)
{
  if (assembly_paths.size() != assembly_indices.size() ||
      assembly_paths.size() != is_targets.size()) {
    throw std::runtime_error(
      "assembly_path, assembly_idx, and is_target must have the same length");
  }

  std::string kmers;
  py::list all_idx_to_id;
  std::vector<std::size_t> record_offsets;
  std::vector<std::size_t> assembly_offsets;
  std::size_t global_minimizer_idx = 0;

  for (std::size_t assembly_i = 0; assembly_i < assembly_paths.size(); ++assembly_i) {
    const auto assembly_idx = assembly_indices[assembly_i];
    if (assembly_idx > std::numeric_limits<uint16_t>::max()) {
      throw std::runtime_error("assembly_idx must fit in uint16");
    }

    const auto records = btllib::read_fasta(assembly_paths[assembly_i]);
    py::tuple idx_to_id(records.size());
    kmers.reserve(kmers.size() + (records.size() * 17));
    bool assembly_has_minimizers = false;

    for (std::size_t record_idx = 0; record_idx < records.size(); ++record_idx) {
      if (record_idx > std::numeric_limits<uint16_t>::max()) {
        throw std::runtime_error("record_idx must fit in uint16");
      }

      const auto& record = records[record_idx];
      idx_to_id[record_idx] = py::str(record.id);

      const auto mins = btllib::minimize_sequence(record.sequence, kmerlen, windowsize);
      if (!mins.empty()) {
        if (!assembly_has_minimizers) {
          assembly_offsets.push_back(global_minimizer_idx);
          assembly_has_minimizers = true;
        }
        record_offsets.push_back(global_minimizer_idx);
      }
      for (const auto& m : mins) {
        if (m.pos > std::numeric_limits<uint32_t>::max()) {
          throw std::runtime_error("minimizer position exceeds uint32 range");
        }

        append_scalar<uint64_t>(kmers, m.out_hash);
        append_scalar<uint32_t>(kmers, static_cast<uint32_t>(m.pos));
        append_scalar<uint16_t>(kmers, static_cast<uint16_t>(record_idx));
        append_scalar<uint16_t>(kmers, static_cast<uint16_t>(assembly_idx));
        append_scalar<uint8_t>(kmers, is_targets[assembly_i] ? uint8_t{1} : uint8_t{0});
        ++global_minimizer_idx;
      }
    }

    all_idx_to_id.append(idx_to_id);
  }

  return { kmers, all_idx_to_id, record_offsets, assembly_offsets };
}

} // namespace

PYBIND11_MODULE(_core, m)
{
  m.doc() = "Minimal btllib indexlr bindings";

  m.def("indexlr_native",
        [](const std::vector<std::string>& assembly_paths,
           std::size_t kmerlen,
           std::size_t windowsize,
           const std::vector<std::size_t>& assembly_indices,
           const std::vector<bool>& is_targets) {
          auto [kmers, ids, record_offsets, assembly_offsets] =
            indexlr_impl(assembly_paths, kmerlen, windowsize, assembly_indices, is_targets);
          return py::make_tuple(
            py::bytes(kmers), ids, py::cast(record_offsets), py::cast(assembly_offsets));
        },
        py::arg("assembly_path"),
        py::arg("kmerlen"),
        py::arg("windowsize"),
        py::arg("assembly_idx"),
        py::arg("is_target"));
}
