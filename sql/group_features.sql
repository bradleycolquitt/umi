SELECT DISTINCT grouped.bc, grouped.hpos1, gene_id, sum(unique_umi), sum(total_umi)
FROM
    grouped JOIN merge ON
        grouped.bc = merge.bc
        AND grouped.tid = merge.tid
        AND grouped.hpos1 = merge.hpos1
        AND grouped.strand = merge. strand
    JOIN maker_anno_run3 ON
        merge.instrument = maker_anno_run3.instrument
        AND merge.flowcell = maker_anno_run3.flowcell
        AND merge.cluster = maker_anno_run3.cluster

GROUP BY
    gene_id

;
