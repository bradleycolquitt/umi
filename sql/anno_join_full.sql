INSERT INTO %s
SELECT
    merge.bc,
    merge.chrom,
    merge.lpos1,
    merge.lpos2,
    merge.rpos1,
    merge.rpos2,
    merge.strand,
    merge.umi,
    -- count(merge.umi) as umi
    -- count(case when merge.strand = anno.genes.strand then merge.umi end)                       as total_umi_same_strand,
    anno.genes.gene_id,
    anno.genes.transcript_id,
    anno.genes.strand,
    anno.genes.element
FROM
   --merge, pos_rtree ON merge.rowid = pos_rtree.id
    anno.genes JOIN merge
        ON merge.chrom = anno.genes.chrom
            AND anno.genes.build=:build
        JOIN pos_rtree ON merge.rowid = pos_rtree.id
   --          AND anno.genes.build=:build
            AND
            (
            pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end OR
                 pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end
            )
--       WHERE anno.genes.build=:build
             -- GROUP BY
    --               merge.bc,
    --               merge.umi,
    -- merge.chrom,
    -- CASE WHEN merge.strand IS 0 THEN merge.lpos1 ELSE merge.rpos2 END
