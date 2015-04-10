               SELECT
                   merge.bc,
                   anno.genes.chrom,
                   merge.hpos,
                   merge.tpos,
                   count(merge.umi) as total_umi,
                   count(case when merge.strand = anno.genes.strand then merge.umi end)
                       as total_umi_same_strand,
                   count(distinct merge.umi) as unique_umi,
                   anno.genes.gene_id,
                   anno.genes.transcript_id,
                   anno.genes.element
               FROM
                   merge, pos_rtree WHERE merge.rowid = pos_rtree.id
                   JOIN reference ON merge.tid = reference.tid
                   JOIN anno.genes ON reference.name = anno.genes.chrom
                        AND anno.genes.build=:build
                        AND
                        (
                        pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end OR
                        pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end
                        )
               GROUP BY
                   merge.bc, merge.tid, merge.rpos1
