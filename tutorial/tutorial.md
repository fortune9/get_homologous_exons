# This is a brief tutorial to identify homologous exons from sequence
# alignment

## Programs to used
> 1. <../homologous_exons_by_gene_alignment.pl>
> 2. <../best_homo_exons.pl>

### Step 1. Get all the overlapped exons from the alignments
> ../homologous_exons_by_gene_alignment.pl -a aln/ \\
> -s struct/dsim2_CDS_exon_string.tsv \\
> -s struct/dyak1_CDS_exon_string.tsv \\
> -s struct/dmel5_CDS_exon_string.tsv --minus-stop -o homo_exons.tsv

### Step 2: find the best homologous exon for each exon in each
sequence
> ../best_homo_exons.pl -o best_exons.tsv -i 70 \\
> --alen-frac 0.9 --bound-dist 10 homo_exons.tsv

## More detailed tutorial is coming.

