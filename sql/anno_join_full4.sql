
SELECT
    merge.bc as bc,
    merge.chrom as chrom,
    count(distinct merge.umi) as unique_umi,
    count(distinct case when collapsed.strand = anno.genes.strand then merge.umi end)
        AS unique_umi_strand,
    anno.genes.gene_id AS gene_id,
    anno.genes.transcript_id AS transcript_id,
    anno.genes.strand AS gene_strand,
    anno.genes.element AS element
FROM
    anno.genes JOIN merge ON merge.chrom = anno.genes.chrom AND anno.genes.build=:build
        JOIN collapsed ON merge.idcollapsed=collapsed.rowid
        AND collapsed.isize BETWEEN 150 AND 1000
        JOIN pos_rtree ON collapsed.rowid = pos_rtree.id
        -- JOIN anno.genes
        --     ON collapsed.chrom = anno.genes.chrom
        --     AND anno.genes.build=:build
            AND
        (
            (pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end)
            OR
            (pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end)
        )

GROUP BY
    merge.strand,
    merge.bc,
    merge.chrom,
    merge.hpos1;
