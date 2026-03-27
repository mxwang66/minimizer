#pragma once

#include <string>
#include <vector>

namespace btllib {

struct FastaRecord
{
  std::string id;
  std::string sequence;
};

std::vector<FastaRecord>
read_fasta(const std::string& assembly_path);

} // namespace btllib
