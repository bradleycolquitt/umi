ATTACH DATABASE '/media/data/db/anno.db' as anno;

SELECT bc, gene_id, sum(unique_umi) as unique_umi, sum(total_umi) as total_umi FROM
-- Identify intersecting read positions
(
SELECT DISTINCT collapsed.bc, anno.genes_ref.gene_id, grouped.unique_umi, grouped.total_umi FROM
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
       JOIN grouped ON collapsed.pos_id=grouped.pos_id
       WHERE
           anno.genes_ref.build='maker'
           --AND anno.genes_ref.tid=0
           --AND anno.genes_ref.tid BETWEEN 0 and 10
) as tmp
GROUP by bc, gene_id;
