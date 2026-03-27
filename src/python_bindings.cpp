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

std::tuple<std::string, py::tuple>
indexlr_impl(const std::string& assembly_path,
             std::size_t kmerlen,
             std::size_t windowsize,
             std::size_t assembly_idx,
             bool is_target)
{
  if (assembly_idx > std::numeric_limits<uint16_t>::max()) {
    throw std::runtime_error("assembly_idx must fit in uint16");
  }

  const auto records = btllib::read_fasta(assembly_path);
  py::tuple idx_to_id(records.size());

  std::string kmers;
  kmers.reserve(records.size() * 17);

  for (std::size_t record_idx = 0; record_idx < records.size(); ++record_idx) {
    if (record_idx > std::numeric_limits<uint16_t>::max()) {
      throw std::runtime_error("record_idx must fit in uint16");
    }

    const auto& record = records[record_idx];
    idx_to_id[record_idx] = py::str(record.id);

    const auto mins = btllib::minimize_sequence(record.sequence, kmerlen, windowsize);
    for (const auto& m : mins) {
      if (m.pos > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("minimizer position exceeds uint32 range");
      }

      append_scalar<uint64_t>(kmers, m.out_hash);
      append_scalar<uint32_t>(kmers, static_cast<uint32_t>(m.pos));
      append_scalar<uint16_t>(kmers, static_cast<uint16_t>(record_idx));
      append_scalar<uint16_t>(kmers, static_cast<uint16_t>(assembly_idx));
      append_scalar<uint8_t>(kmers, is_target ? uint8_t{1} : uint8_t{0});
    }
  }

  return { kmers, idx_to_id };
}

} // namespace

PYBIND11_MODULE(_core, m)
{
  m.doc() = "Minimal btllib indexlr bindings";

  m.def("indexlr_native",
        [](const std::string& assembly_path,
           std::size_t kmerlen,
           std::size_t windowsize,
           std::size_t assembly_idx,
           bool is_target) {
          auto [kmers, ids] =
            indexlr_impl(assembly_path, kmerlen, windowsize, assembly_idx, is_target);
          return py::make_tuple(py::bytes(kmers), ids);
        },
        py::arg("assembly_path"),
        py::arg("kmerlen"),
        py::arg("windowsize"),
        py::arg("assembly_idx"),
        py::arg("is_target"));
}
