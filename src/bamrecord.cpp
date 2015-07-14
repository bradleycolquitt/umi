#include <bamrecord.h>

BamRecord::BamRecord(int tid)
    : read_id("")
    , gene_id("")
    , tid(tid)
    , pos(0)
    , bc(-1)
    , umi(0)
    , umi2("AA")
{
}

BamRecord::BamRecord(const string read_id, const string gene_id)
    : read_id(read_id)
    , gene_id(gene_id)
    , tid(0)
    , pos(-1)
    , bc(-1)
    , umi(0)
    , umi2("AA")
{
}

void BamRecord::set_bc(BamHash* bamhash, bam1_t* b, int* used_offset) {
    bc = get_sequence(b, bamhash->get_sequence_pos(2), bamhash->get_sequence_pos(3), bamhash->get_barcodes(), bamhash->get_bc_offsets(), bamhash->get_bc_min_qual(), used_offset);
}

void BamRecord::set_umi(BamHash* bamhash, bam1_t* b, int used_offset) {
    umi = get_sequence(b, bamhash->get_sequence_pos(0), bamhash->get_sequence_pos(1), used_offset);
}

void BamRecord::set_bc(BamHash* bamhash, char* seq, char* qual, int* used_offset) {
    bc = get_sequence(seq, qual, bamhash->get_sequence_pos(2), bamhash->get_sequence_pos(3), bamhash->get_barcodesA(), bamhash->get_bc_offsets(), bamhash->get_bc_min_qual(), used_offset);
}

void BamRecord::set_umi(BamHash* bamhash, char* seq, char* qual, int used_offset) {
    umi2 = get_sequence(seq, qual, bamhash->get_sequence_pos(0), bamhash->get_sequence_pos(1), bamhash->get_bc_min_qual(), used_offset);
}

bool BamRecord::is_complete() {
    return (pos>0 && bc>-1);
}
