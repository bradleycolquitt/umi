INSERT INTO %s
SELECT
    collapsed.bc,
    collapsed.chrom,
    collapsed.lpos1,
    collapsed.lpos2,
    collapsed.rpos1,
    collapsed.rpos2,
    collapsed.strand,
    collapsed.total_umi,
    --count(collapsed.umi) as total_umi
    --count(case when collapsed.strand = anno.genes.strand then collapsed.umi end)
    anno.genes.gene_id,
    anno.genes.transcript_id,
    anno.genes.strand,
    anno.genes.element
FROM
    collapsed, pos_rtree ON collapsed.rowid = pos_rtree.id
    JOIN anno.genes
        ON collapsed.chrom = anno.genes.chrom
        AND anno.genes.build=:build
        AND
        (
            pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end OR
            pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end
        )
