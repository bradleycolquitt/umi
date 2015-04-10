ATTACH DATABASE "merge.db" AS merge_db;

CREATE TABLE IF NOT EXISTS merge_db.merge (instrument text, flowcell text, cluster text, chrom text, tid int, lpos1 int, lpos2 int, rpos1, rpos2, strand int, bc int, umi int);

INSERT INTO merge_db.merge
    SELECT
        read1.instrument,
        read1.flowcell,
        read1.cluster,
        read1.chrom,
        read1.tid,
        CASE WHEN read1.strand IS 0 THEN read1.hpos ELSE read2.hpos END,
        CASE WHEN read1.strand IS 0 THEN read1.tpos ELSE read2.tpos END,
        CASE WHEN read1.strand IS 0 THEN read2.tpos ELSE read1.tpos END,
        CASE WHEN read1.strand IS 0 THEN read2.hpos ELSE read1.hpos END,
        read1.strand,
        read1.bc,
        read1.umi
    FROM read1 LEFT JOIN read2 ON read1.cluster=read2.cluster;
