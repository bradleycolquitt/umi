#include <bam_utils.h>
#include <bam_sql_merge.h>
#include <htslib/sam.h>
using namespace std;

int main(int argc, char ** argv) {

   const char * bam_fname = "/media/data/bam/150403/lonStrDom1/lib40/lib40_clip18Aligned.bam";
   samFile* bam = sam_open(bam_fname, "rb");
   bam_hdr_t* header = sam_hdr_read(bam);
   hts_idx_t* idx = bam_index_load(bam_fname);

   hts_itr_t* bam_itr;
   bam_itr = bam_itr_queryi(idx, 0, 0, header->target_len[0]);
   bam1_t* b = bam_init1();

   BamDB* bamdb = new BamDB(bam_fname, "dummy.db", "bc8", 8, 17, true);
   int i = 0;

   map<int,int> bc_counts;
   for (int i = -2; i < 8; ++i) {
       bc_counts[i] = 0;
   }

   int result = 0;
   while ((result = sam_itr_next(bam, bam_itr, b)) >= 0 ) {
       int seqlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
       if (seqlen > 100) continue;
       //if (i > 5) break;
       if (bam_is_rev(b) && ((b->core.flag&BAM_FREAD1) != 0)) {
           ++i;
           //cout << b->core.flag << endl;
           //cout << (b->core.flag&BAM_FREAD1) << endl;
           uint8_t* seq = bam_get_seq(b);
           uint8_t* qual = bam_get_qual(b);
           uint8_t seq2[seqlen];
           uint8_t qual2[seqlen];

           bam_get_seq_qual(b, seq2, qual2);

           //print_uint8(seq, seqlen, true);
           //print_uint8(seq2, seqlen, false);
           //cout << endl;
           int used_offset = 0;
           int bc = get_sequence(b, 8, 13, bamdb->get_barcodes(), bamdb->get_bc_offsets(), 10, &used_offset);
           //cout << bc << endl;
           bc_counts[bc] = bc_counts[bc] + 1;
       }

    }
    for (int i = -2; i < 8 ; ++i) {
       cout << "bc" << i <<": " << bc_counts[i] << endl;
   }
}
