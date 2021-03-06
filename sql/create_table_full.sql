CREATE TABLE IF NOT EXISTS %s
                   (bc int,
                   chrom int,
                   lpos1 int,
                   lpos2 int,
                   rpos1 int,
                   rpos2 int,
                   strand int,
                   umi int,
                   -- umi_count int,
                   -- total_umi_same_strand int,
                   gene_id text,
                   transcript_id text,
                   gene_strand text,
                   element text);
