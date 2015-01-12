### Log
#### 1/9/15
Comparison between raw featureCounts and either sum of filtered (bc and umi quality) or binomial MLE of n indicates that raw counts are better than fitlered or inferred expression estimate. 

So will hold off on development of read db for expression quant
However, bam_sql_cy2 will still be potentially useful for quick access to expression levels.

Look into just bc filtering compared to raw counts.

Most filtering due to cigar filter

#### 1/11/15
Major descrepancy between featureCounts and load_to_db2 is use of full read in featureCounts. This produces much better correlation with ERCC molecule number. Using --read2pos 5 flag of featureCounts recapitulates pattern seen in load_to_db2 of overestimatation of expression of lowly adundant transcripts. Perhaps this is due to slight alignment errors and is simply an artifcat of having the transcripts concatenated.

Will make ERCC genome with 100 bp spacers. DONE. saved as ERCC92_cat_space.fa with corresponding .gtf in lib/gtf.

Filtering out reads with tags containing [IDN] fixes problem! Removes split reads. Read counts very well correlated with transcript number now (R.squared ~0.96).

Might be issue with reads spanning introns in actual data. Will have to look out for.

