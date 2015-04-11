#include <dbrecord.h>
//#include <bam_utils.h>
//#include <bam_sql_merge.h>
//#include <string_utils.h>
//#include <sqlite3.h>
//#include <cstring>

dbRecord::dbRecord(BamDB* bamdb)
    : bamdb(bamdb)
{
    prepare_exists();
}

dbRecord::~dbRecord() {}

void dbRecord::prepare_exists() {
    int result = 0;
    const char* tail = 0;
    char SQL[BUFFER_SIZE];
    sprintf(SQL, "SELECT EXISTS(SELECT 1 FROM align WHERE cluster=?);");
    if ((result = sqlite3_prepare_v2(bamdb->get_conn(), SQL, BUFFER_SIZE, \
                                     &stmt_exists, &tail)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on record check prep. \n", result);
    }
}


/***********************************
Read processing and record updating
************************************/
// Split qname/read metadata
void dbRecord::split_qname(bam1_t* b) {
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

void dbRecord::set_bc(BamDB* bamdb, bam1_t* b, int* used_offset) {
    bc = get_sequence(b, bamdb->get_sequence_pos(2), bamdb->get_sequence_pos(3), bamdb->get_barcodes(), bamdb->get_bc_offsets(), bamdb->get_bc_min_qual(), used_offset);
}

void dbRecord::set_umi(BamDB* bamdb, bam1_t* b, int used_offset) {
    umi = get_sequence(b, bamdb->get_sequence_pos(0), bamdb->get_sequence_pos(1), used_offset);
}

dbRecordSe::dbRecordSe(BamDB* bamdb)
    : dbRecord(bamdb)
{
    read_pos.resize(2);
    //Prepare insert statement
    const char* tail = 0;
    char SQL[BUFFER_SIZE];
    sprintf(SQL, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS, @TPOS, @STR, @BC, @UMI);");
    sqlite3_prepare_v2(bamdb->get_conn(), SQL, BUFFER_SIZE, &stmt, &tail);
}
int dbRecordSe::set_positions(bam1_t* b, int read_num) {
            if (bam_is_rev(b)) {
                read_pos[1] = b->core.pos + 1;
                read_pos[0] = bam_endpos(b);
                strand = true;
            } else {
                read_pos[0] = b->core.pos + 1;
                read_pos[1] = bam_endpos(b);
                strand = false;
            }
            return 0;
}

int dbRecordSe::insert_to_db() {
    sqlite3_bind_text(stmt, 1, instrument, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, flowcell, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, cluster, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(stmt, 4, tid);
    sqlite3_bind_int(stmt, 5, read_pos[0]);
    sqlite3_bind_int(stmt, 6, read_pos[1]);
    sqlite3_bind_int(stmt, 7, strand);
    sqlite3_bind_int(stmt, 8, bc);
    sqlite3_bind_int(stmt, 9, umi);

    sqlite3_step(stmt);
    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);

    return 0;
}

dbRecordPe::dbRecordPe(BamDB* bamdb)
: dbRecord(bamdb)
{
     read_pos.resize(4);
    int result = 0;
    const char* tail = 0;

    char insert_sql[BUFFER_SIZE];
    sprintf(insert_sql, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS1, @TPOS1, @HPOS2, @TPOS2, @INS, @STR, @BC, @UMI);");
    sqlite3_prepare_v2(bamdb->get_conn(), insert_sql, BUFFER_SIZE, \
                       &insert_stmt, &tail);

    //Prepare update statement, read1
    char update_sql1[BUFFER_SIZE];
    sprintf(update_sql1, "UPDATE align SET hpos1=?, tpos1=?, strand=?, bc=?, umi=? WHERE cluster=?");
    sqlite3_prepare_v2(bamdb->get_conn(), update_sql1, BUFFER_SIZE, \
                       &update_stmt1, &tail);

    //Prepare update statement, read2
    char update_sql2[BUFFER_SIZE];
    sprintf(update_sql2, "UPDATE align SET hpos2=?, tpos2=?, isize=? WHERE cluster=?;");
    if ((result = sqlite3_prepare_v2(bamdb->get_conn(), update_sql2, BUFFER_SIZE,\
                                     &update_stmt2, &tail)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on update2 prep. \n", result);
    }

}

int dbRecordPe::set_positions(bam1_t* b, int read_num) {
    if (read_num == 1) {
    if (bam_is_rev(b)) {
            read_pos[0] = bam_endpos(b);
            read_pos[1] = b->core.pos + 1;
            strand = true;
        } else {
            read_pos[0] = b->core.pos + 1;
            read_pos[1] = bam_endpos(b);
            strand = false;
        }
    } else if (read_num == 2) {
        if (bam_is_rev(b)) {
            read_pos[2] = bam_endpos(b);
            read_pos[3] = b->core.pos + 1;
        } else {
            read_pos[2] = b->core.pos + 1;
            read_pos[3] = bam_endpos(b);
        }
    }
    return 0;
}

int dbRecordPe::exists(BamDB* bamdb) {
    int result = 0;
    int exec_result = 0;

    // #ifdef DEBUG
    // cerr << cluster << endl;
    // #endif
    DEBUG_LOG("cluster:" << cluster);
    sqlite3_bind_text(stmt_exists, 1, cluster, -1, SQLITE_TRANSIENT);
    if ((result = sqlite3_step(stmt_exists)) < 100) {
        fprintf(stderr, "SQL error (%d) on record check execution. \n", result);
    }
    exec_result = sqlite3_column_int(stmt_exists, 0);
    DEBUG_LOG("exists result:" << exec_result);
    sqlite3_reset(stmt_exists);
    return(exec_result);
}

int dbRecordPe::insert_to_db() {
    return 0;
}

int dbRecordPe::insert_to_db(int read_num) {
    int result = 0;
    // #ifdef DEBUG
    //      cerr << "insert read" << read_num << endl;
    // //     cerr << "hpos1:" << read_pos[0]
    // //          << " tpos1:" << read_pos[1] << endl;
    // #endif
    //DEBUG_LOG("y");

    if (read_num == 1) {
        sqlite3_bind_text(insert_stmt, 1, instrument, -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 2, flowcell, -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 3, cluster, -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(insert_stmt, 4, tid);
        sqlite3_bind_int(insert_stmt, 5, read_pos[0]);
        sqlite3_bind_int(insert_stmt, 6, read_pos[1]);
        sqlite3_bind_int(insert_stmt, 7, 0);
        sqlite3_bind_int(insert_stmt, 8, 0);
        sqlite3_bind_int(insert_stmt, 9, 0);
        sqlite3_bind_int(insert_stmt, 10, strand);
        sqlite3_bind_int(insert_stmt, 11, bc);
        sqlite3_bind_int(insert_stmt, 12, umi);
    } else if (read_num == 2) {
        sqlite3_bind_text(insert_stmt, 1, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 2, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 3, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(insert_stmt, 4, 0);
        sqlite3_bind_int(insert_stmt, 5, 0);
        sqlite3_bind_int(insert_stmt, 6, 0);
        sqlite3_bind_int(insert_stmt, 7, read_pos[2]);
        sqlite3_bind_int(insert_stmt, 8, read_pos[3]);
        sqlite3_bind_int(insert_stmt, 9, insert);
        sqlite3_bind_int(insert_stmt, 10, 0);
        sqlite3_bind_int(insert_stmt, 11, 0);
        sqlite3_bind_int(insert_stmt, 12, 0);
    }
    if ((result = sqlite3_step(insert_stmt)) != SQLITE_DONE ) {
        fprintf(stderr, "Insertion error (%d): read_num=%d, %s\n", result, \
               read_num, cluster);
    }
    sqlite3_clear_bindings(insert_stmt);
    sqlite3_reset(insert_stmt);
    return 0;
}

int dbRecordPe::update_record(int read_num) {
    int result = 0;
    // #ifdef DEBUG
    //     cerr << "update read" << read_num << endl;
    // #endif
    if (read_num == 1) {
        sqlite3_bind_int(update_stmt1, 1, read_pos[0]);
        sqlite3_bind_int(update_stmt1, 2, read_pos[1]);
        sqlite3_bind_int(update_stmt1, 3, strand);
        sqlite3_bind_int(update_stmt1, 4, bc);
        sqlite3_bind_int(update_stmt1, 5, umi);
        sqlite3_bind_text(update_stmt1, 6, cluster, -1, SQLITE_TRANSIENT);

        if ((result = sqlite3_step(update_stmt1)) != SQLITE_DONE ) {
            fprintf(stderr, "Update error (%d): read_num=1, %s\n", result, cluster);
        }
        sqlite3_clear_bindings(update_stmt1);
        sqlite3_reset(update_stmt1);
    } else if (read_num == 2) {
        sqlite3_bind_int(update_stmt2, 1, read_pos[2]);
        sqlite3_bind_int(update_stmt2, 2, read_pos[3]);
        sqlite3_bind_int(update_stmt2, 3, insert);
        sqlite3_bind_text(update_stmt2, 4, cluster, -1, SQLITE_TRANSIENT);

        if ((result = sqlite3_step(update_stmt2)) != SQLITE_DONE ) {
            fprintf(stderr, "Update error (%d): read_num=2, %s\n", result, cluster);
        }
        sqlite3_clear_bindings(update_stmt2);
        sqlite3_reset(update_stmt2);
    }
    return 0;
}

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
    // DEBUG_LOG(instrument);
    // DEBUG_LOG(flowcell);
    // DEBUG_LOG(chrom);
    // DEBUG_LOG(tid);
    // DEBUG_LOG(strand);
    // DEBUG_LOG(read_pos[0]);
    // DEBUG_LOG(read_pos[1]);
    // DEBUG_LOG(bc);
    // DEBUG_LOG(umi);


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
    //if ((result = sqlite3_step(insert_stmt)) != SQLITE_DONE ) {
    if (result != SQLITE_DONE) {
        throw sql_exception(result, sqlite3_errmsg(bamdb->get_conn()));
        //fprintf(stderr, "Insertion error (%d): read1, %s\n", result, cluster);

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
    sprintf(insert_sql, "INSERT INTO read2 VALUES (@IN, @FL, @CL, @CHR, @TID, @HPOS, @TPOS, @STR);");
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

    if ((result = sqlite3_step(insert_stmt)) != SQLITE_DONE ) {
        throw sql_exception(result, sqlite3_errmsg(bamdb->get_conn()));
        //fprintf(stderr, "Insertion error (%d): read2, %s\n", result, cluster);
    }
    sqlite3_clear_bindings(insert_stmt);
    sqlite3_reset(insert_stmt);
    return 0;
}
