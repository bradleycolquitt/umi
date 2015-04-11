INSERT INTO %s
SELECT
    collapsed.bc as bc,
    collapsed.chrom as chrom,
    count(distinct merge.umi) as unique_umi,
    count(distinct case when collapsed.strand = anno.genes.strand then merge.umi end)
        AS unique_umi_strand,
    anno.genes.gene_id AS gene_id,
    anno.genes.transcript_id AS transcript_id,
    anno.genes.strand AS gene_strand,
    anno.genes.element AS element
FROM
    merge JOIN collapsed ON merge.idcollapsed=collapsed.rowid
        JOIN pos_rtree ON collapsed.rowid = pos_rtree.id
        JOIN anno.genes
            ON collapsed.chrom = anno.genes.chrom
            AND anno.genes.build=:build
            AND
        (
            (pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end) OR
            (pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end)
        )
GROUP BY
    collapsed.bc,
    collapsed.chrom,
    CASE WHEN collapsed.strand IS 0 THEN merge.lpos1 ELSE merge.rpos2 END,
    collapsed.strand
