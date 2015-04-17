#include <bamrecord.h>

BamRecord::BamRecord(int tid)
    : read_id("")
    , gene_id("")
    , tid(tid)
    , pos(0)
    , bc(-1)
    , umi(0)
{
}


void BamRecord::set_bc(BamHash* bamhash, bam1_t* b, int* used_offset) {
    bc = get_sequence(b, bamhash->get_sequence_pos(2), bamhash->get_sequence_pos(3), bamhash->get_barcodes(), bamhash->get_bc_offsets(), bamhash->get_bc_min_qual(), used_offset);
}

void BamRecord::set_umi(BamHash* bamhash, bam1_t* b, int used_offset) {
    umi = get_sequence(b, bamhash->get_sequence_pos(0), bamhash->get_sequence_pos(1), used_offset);
}
