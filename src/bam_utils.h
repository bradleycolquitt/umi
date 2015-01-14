#ifndef BAM_UTIL_H
#define BAM_UTIL_H

uint64_t count_bam_records(hts_idx_t* idx, bam_hdr_t* header);

vector<int> seq2int(string& s);

#endif
