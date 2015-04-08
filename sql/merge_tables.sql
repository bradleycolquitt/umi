CREATE TABLE IF NOT EXISTS merge (instrument text, flowcell text, cluster text, tid int, hpos1 int, tpos1 int, hpos2, tpos2, strand int, bc int, umi int);

INSERT INTO merge
    SELECT
        read1.instrument,
        read1.flowcell,
        read1.cluster,
        read1.tid,
        read1.hpos,
        read1.tpos,
        read2.hpos,
        read2.tpos,
        read1.strand,
        read1.bc,
        read1.umi
    FROM read1 LEFT JOIN read2 ON read1.cluster=read2.cluster;
