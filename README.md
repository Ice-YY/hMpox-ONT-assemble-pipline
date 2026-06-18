# hMpox-ONT-assemble-pipline
A pipeline used by Hangzhou CDC for assembling ONT raw data against a reference sequence.
## Requirements

- Python 3
- `fastplong`, `minimap2`, `samtools`, `bcftools`
- Docker (for Clair3)

## Usage

python run_pipeline.py --raw sample.fastq.gz --ref ref.fasta

## Reference 

1. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890. doi:10.1093/bioinformatics/bty560
2. Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191
3. Zheng Z, Li S, Su J, Leung AW, Lam T, Luo R. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. Nat Comput Sci. 2022;2(12):797-803. doi:10.1038/s43588-022-00387-x
4. Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008
