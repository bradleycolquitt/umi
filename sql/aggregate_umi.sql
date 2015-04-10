SELECT
    bc,
    chrom,
    lpos1,
    lpos2,
    rpos1,
    rpos2,
    umi,
    count(umi),
    count(case when strand = gene_strand then umi end),
    gene_id,
    transcript_id,
    gene_strand,
    element
FROM
    %s
GROUP BY
    bc,
    umi,
    chrom,
    CASE WHEN strand IS 0 THEN lpos1 ELSE rpos2 END
