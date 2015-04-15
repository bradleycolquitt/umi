ATTACH DATABASE '/media/data/db/anno.db' as anno;

--SELECT tmp.bc as bc, tmp.gene_id as gene_id, sum(tmp.unique_umi) as unique_umi, sum(tmp.total_umi) as total_umi FROM
--(
--SELECT bc, gene_id, sum(unique_umi) as unique_umi, sum(total_umi) as total_umi FROM

--SELECT * FROM
CREATE TEMP TABLE tmp AS SELECT distinct collapsed.bc as bc, anno_tid.strand, anno_tid.gene_id as gene_id, grouped.pos_id, grouped.unique_umi as unique_umi, grouped.total_umi as total_umi FROM
--SELECT collapsed.bc as bc, collapsed.anno_tid.gene_id as gene_id, sum(grouped.unique_umi) as unique_umi, sum(grouped.total_umi) as total_umi FROM

-- Associate chrom with tid
(SELECT * FROM anno.genes JOIN reference ON anno.genes.chrom=reference.name and anno.genes.build='maker') AS anno_tid

-- Identify intersecting read positions
LEFT OUTER JOIN
  pos_rtree ON
      anno_tid.tid=pos_rtree.tid1
      AND pos_rtree.strand2=anno_tid.strand
      AND (
          pos_rtree.lpos1 BETWEEN anno_tid.start AND anno_tid.end
          OR pos_rtree.rpos2 BETWEEN anno_tid.start AND anno_tid.end
          )

-- Join to grouped
JOIN collapsed ON pos_rtree.id=collapsed.rowid
JOIN grouped ON collapsed.pos_id=grouped.pos_id;

CREATE INDEX tmp_idx ON tmp(bc, gene_id);
SELECT bc, gene_id, sum(unique_umi) as unique_umi, sum(total_umi) as total_umi FROM tmp GROUP BY bc, gene_id
--LIMIT 0,6;
-- Group by pos_id, count number of distinct hpos2 (rpos2/lpos1) associated with pos_id, filter on this value
--GROUP BY grouped.pos_id
--HAVING count(*) > 10
--)
--GROUP BY bc, gene_id
-- GROUP BY collapsed.bc, anno_tid.gene_id
-- HAVING count(grouped.pos_id) > 2
;
