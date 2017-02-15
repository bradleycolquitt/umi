SELECT
    tid,
    hpos1,
    strand,
    count(distinct umi),
    count(umi),
FROM
    merge
GROUP BY
    tid,
    hpos1,
    strand
