DROP TABLE IF EXISTS collapsed;

CREATE TABLE IF NOT EXISTS collapsed (
    bc int,
    umi int,
    chrom text,
    lpos1 int,
    lpos2 int,
    rpos1 int,
    rpos2 int,
    strand int,
    total_umi int
);

INSERT INTO collapsed
SELECT
    bc,
    umi,
    chrom,
    lpos1,
    lpos2,
    rpos1,
    rpos2,
    strand,
    count(umi)
FROM
    merge
GROUP BY
    bc,
    umi,
    chrom,
    strand,
    CASE WHEN strand IS 0 THEN lpos1 ELSE rpos2 END
;
