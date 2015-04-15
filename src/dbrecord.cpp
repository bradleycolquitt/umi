#include <dbrecord.h>
//#include <bam_utils.h>
//#include <bam_sql_merge.h>
//#include <string_utils.h>
//#include <sqlite3.h>
//#include <cstring>



void dbRecord0::split_qname(bam1_t* b) {
    char* qname = bam_get_qname(b);
    char* qname_array = strtok(qname, ":");
    int j = 0;
    while (qname_array != NULL) {
        if (j == 0) {
            strcpy(instrument, qname_array);
        } else if (j == 2) {
            strcpy(flowcell, qname_array);
        } else if (j== 3) {
            strcpy(cluster, qname_array);
        } else if (j > 3) {
            strcat(cluster, ":");
            strcat(cluster, qname_array);
        }
        qname_array = strtok(NULL, ":");
        ++j;
    }
}

int dbRecord0::set_positions(bam1_t* b) {
    if (bam_is_rev(b)) {
        read_pos[0] = bam_endpos(b);
        read_pos[1] = b->core.pos + 1;
        strand = true;
    } else {
        read_pos[0] = b->core.pos + 1;
        read_pos[1] = bam_endpos(b);
        strand = false;
    }
}

void dbRecord0::set_tid(int32_t tid) {
    this->tid = tid;
}

void dbRecord0::set_chrom(BamDB* bamdb, int32_t tid){
    strcpy(chrom, bamdb->get_chrom(tid));
}

dbRecord1::dbRecord1(BamDB* bamdb)
: dbRecord0(bamdb)
{
    const char* tail = 0;

    char insert_sql[BUFFER_SIZE];
    sprintf(insert_sql, "INSERT INTO read1 VALUES (@IN, @FL, @CL, @CHR, @TID, @HPOS, @TPOS, @STR, @BC, @UMI);");
    sqlite3_prepare_v2(bamdb->get_conn(), insert_sql, BUFFER_SIZE, \
                       &insert_stmt, &tail);
}

void dbRecord1::set_bc(BamDB* bamdb, bam1_t* b, int* used_offset) {
    bc = get_sequence(b, bamdb->get_sequence_pos(2), bamdb->get_sequence_pos(3), bamdb->get_barcodes(), bamdb->get_bc_offsets(), bamdb->get_bc_min_qual(), used_offset);
}

void dbRecord1::set_umi(BamDB* bamdb, bam1_t* b, int used_offset) {
    umi = get_sequence(b, bamdb->get_sequence_pos(0), bamdb->get_sequence_pos(1), used_offset);
}

int dbRecord1::insert_to_db() {
    int result = 0;
    sqlite3_bind_text(insert_stmt, 1, instrument, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 2, flowcell, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 3, cluster, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 4, chrom, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 5, tid);
    sqlite3_bind_int(insert_stmt, 6, read_pos[0]);
    sqlite3_bind_int(insert_stmt, 7, read_pos[1]);
    sqlite3_bind_int(insert_stmt, 8, strand);
    sqlite3_bind_int(insert_stmt, 9, bc);
    sqlite3_bind_int(insert_stmt, 10, umi);

    result = sqlite3_step(insert_stmt);
    if (result != SQLITE_DONE) {
        throw sql_exception(result, sqlite3_errmsg(bamdb->get_conn()));
    }

    sqlite3_reset(insert_stmt);
    sqlite3_clear_bindings(insert_stmt);

    return 0;
}

dbRecord2::dbRecord2(BamDB* bamdb)
: dbRecord0(bamdb)
{
    const char* tail = 0;
    char insert_sql[BUFFER_SIZE];
    sprintf(insert_sql, "INSERT INTO read2 VALUES (@IN, @FL, @CL, @CHR, @TID, @HPOS, @TPOS, @STR, @INS);");
    sqlite3_prepare_v2(bamdb->get_conn(), insert_sql, BUFFER_SIZE, \
                       &insert_stmt, &tail);
}

int dbRecord2::insert_to_db() {
    int result = 0;
    sqlite3_bind_text(insert_stmt, 1, instrument, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 2, flowcell, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 3, cluster, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(insert_stmt, 4, chrom, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 5, tid);
    sqlite3_bind_int(insert_stmt, 6, read_pos[0]);
    sqlite3_bind_int(insert_stmt, 7, read_pos[1]);
    sqlite3_bind_int(insert_stmt, 8, strand);
    sqlite3_bind_int(insert_stmt, 9, isize);
    //DEBUG_LOG(isize);
    if ((result = sqlite3_step(insert_stmt)) != SQLITE_DONE ) {
        throw sql_exception(result, sqlite3_errmsg(bamdb->get_conn()));
    }
    sqlite3_clear_bindings(insert_stmt);
    sqlite3_reset(insert_stmt);
    return 0;
}

void dbRecord2::set_isize(int32_t is) {
    this->isize = is;
}
