ATTACH DATABASE '/media/data/db/anno.db' as anno;

--SELECT tmp.bc as bc, tmp.gene_id as gene_id, sum(tmp.unique_umi) as unique_umi, sum(tmp.total_umi) as total_umi FROM
--(
SELECT bc, lpos1, rpos2, gene_id, sum(unique_umi) as unique_umi, sum(total_umi) as total_umi FROM

-- (
--SELECT DISTINCT collapsed.bc as bc, collapsed.anno_tid.gene_id as gene_id, sum(grouped.unique_umi) as unique_umi, sum(grouped.total_umi) as total_umi FROM

-- Associate chrom with tid
--(
--SELECT * FROM
--(
-- SELECT reference.tid, anno.genes.start, anno.genes.end, anno.genes.strand, anno.genes.gene_id, anno.genes.gene_name
--     FROM anno.genes JOIN reference ON anno.genes.chrom=reference.name and anno.genes.build='maker'
--) AS anno_tid

-- Identify intersecting read positions
--LEFT OUTER JOIN
(
SELECT DISTINCT * FROM
--SELECT anno.genes_ref.start, anno.genes_ref.end, anno.gene_rtree.start, anno.gene_rtree.end, pos_rtree.lpos1, pos_rtree.rpos2 FROM
  anno.genes_ref JOIN
      anno.gene_rtree ON anno.genes_ref.rowid = anno.gene_rtree.id
      LEFT JOIN
      collapsed ON
          anno.genes_ref.tid=collapsed.tid
          AND anno.genes_ref.strand=collapsed.strand
      JOIN pos_rtree
          ON collapsed.rowid = pos_rtree.id
          AND
          (
               pos_rtree.lpos1 BETWEEN anno.gene_rtree.start AND anno.gene_rtree.end
               OR pos_rtree.rpos2 BETWEEN anno.gene_rtree.start AND anno.gene_rtree.end
          )

-- Join to grouped
-- JOIN collapsed ON pos_rtree.id=collapsed.rowid
       JOIN grouped ON collapsed.pos_id=grouped.pos_id
       WHERE anno.genes_ref.tid=100 AND anno.genes_ref.build='maker'
) as tmp
GROUP by bc, gene_id;
--LIMIT 0,6;
-- Group by pos_id, count number of distinct hpos2 (rpos2/lpos1) associated with pos_id, filter on this value
--GROUP BY grouped.pos_id
--HAVING count(*) > 10
--)
--GROUP BY bc, gene_id
-- GROUP BY collapsed.bc, anno_tid.gene_id
-- HAVING count(grouped.pos_id) > 2
;
