#include <htslib/sam.h>
#include <bam_utils.h>
#include <sqlite3>

using namespace std;
// Fill class
class BamDB {

    // Variables
    const char* bam_fname;
    samFile* bam;
    bam_hdr_t* header;
    hts_idx_t* idx;
    unsigned int total_mapped;
    const char* db_fname;
    ppdb* conn;
    vector<string> barcodes;


    BamDB(string bam_fname, string dest_fname, string barcodes)
        : bam_fname(bam_fname)
        , dest_fname(dest_fname)
        {
            bam = sam_open(bam_fname, "rb");
            header = sam_hdr_read(bam);
            idx = bam_index_load(bam_fname);
            total_mapped = count_bam_records(idx, header);



        }
    ~BamDB() {
        free(mapped);
        free(unmapped);
    }


    //fill_db
    vector<string>* get_barcodes() {
        return barcodes;
    }
};

int create_table(BamDB& bamdb) {
    const char* statement = "CREATE TABLE IF NOT EXISTS align (name text,
                                                           instrument text,
                                                           flowcell text,
                                                           chrom int,
                                                           position int,
                                                           strand int,
                                                           bc int,
                                                           umi text);";
    sqlite3_exec(bamdb->conn, statement);
    return 0;
}

int fill_db(BamDB& bamdb) {
    boost::regex pattern = "[IDN]";
    int tid = 0;
    hts_itr_t* bam_itr = bam_itr_queryi(bamdb->idx, tid, 0, bamdb->refs[tid]);
    bam1_t r;

    uint32_t* cigar;
    uint8_t* aux_p;
    char* aux;

    string qname;
    int use = 1;

    uint32_t pos;
    while (bam_itr_next(bamdb->bam->fp.bgzf, bam_itr, r)) {

        // get cigar
        cigar = bam_get_cigar(r);

        for (int i=0; i<r->core.n_cigar; ++i) {
            uint32_t op = bam_cigar_op(cigar[i]);
            if (op > 0 & op < 4) {
                use = 0;
                break;
            }
        }
        if (!use) continue;

        // get aux. filter out NH>1
        aux_p = bam_get_aux(r);
        aux = bam_aux2A(aux_p);
        size_t aux_len = strlen(aux);

        boost::regex pattern = "(NH\\:[0-9]+)";
        boost::match_results<string::const_iterator> what;
        while(boost::regex_search(aux, what, pattern, flags)) {
            char* match = what[0];
            if (strlen(match) > 4) {
                use = 0;
                break;
            }
            int mm = atoi(match[3]);
            if (mm > 1) {
                use = 0;
                break;
            }
        }
        if (!use) continue;

        qname = *bam_get_qname(r);
        vector<string> qname_list = split(qname); // need to include from R project
        instrument = qname_list[0];
        flowcell = qname_list[2];

        if (bam_is_rev(r)) {
            pos = bam_endpos(r);
        } else {
            pos = bam->core.pos;
        }

        seq_p = bam_get_seq(r);
        qual_p = bam_get_qual(r);

        bc = get_barcode(seq_p, qual_p, 5, 11, bamdb->barcodes);
        umi = get_umi(seq_p, qual_p, 0, 4);
    }

    return 0;
}

int get_barcode(uint32_t)
