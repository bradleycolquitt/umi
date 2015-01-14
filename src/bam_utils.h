#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include <vector>
#include <string>
#include <htslib/sam.h>

uint64_t count_bam_records(hts_idx_t* idx, bam_hdr_t* header);

std::vector<int> seq2int(std::string& s);

#endif
