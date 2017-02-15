CREATE INDEX IF NOT EXISTS counts_id ON counts(gene_id,bc,i7);

CREATE VIEW IF NOT EXISTS maxes AS
SELECT gene_id, bc, i7, max(count) AS max_count
FROM counts
GROUP BY gene_id,bc,i7;

DROP VIEW IF EXISTS counts_filtered;

CREATE VIEW counts_filtered as SELECT i5,counts.i7,counts.bc,counts.gene_id,umi,count FROM counts
JOIN maxes WHERE
    counts.gene_id = maxes.gene_id AND
    counts.bc = maxes.bc AND
    counts.i7 = maxes.i7 AND
    counts.count > (maxes.max_count * 0.05);

DROP TABLE IF EXISTS collapse;

CREATE TABLE collapse as SELECT
    i5,
    i7,
    gene_id,
    bc,
    count(*) AS unique_umi,
    sum(count) AS total_umi,
    sum(count) / count(*)
FROM counts_filtered
GROUP BY gene_id,bc,i7;
