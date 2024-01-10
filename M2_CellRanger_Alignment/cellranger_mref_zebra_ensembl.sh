

# Download genome assembly
wget https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz

# Download annotation file
wget https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz
gunzip Danio_rerio.GRCz11.110.gtf.gz

# Filter gtf index from non-polyA transcripts (multi-mapped)
cellranger mkgtf \
  Danio_rerio.GRCz11.110.gtf \
  Danio_rerio.GRCz11.110.ensembl.filtered.gtf \
  --attribute=gene_biotype:protein_coding

# make index genome reference file for cell ranger
cellranger mkref \
  --genome=Danio-rerio-genome-ensembl \
  --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa \
  --genes=Danio_rerio.GRCz11.110.ensembl.filtered.gtf \
  --ref-version="fa-GRCz11_gtf-GRCz11.110-ensembl" \
  --nthreads=12
