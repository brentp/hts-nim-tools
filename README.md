# hts-nim-tools

This repository contains a number of tools created with [hts-nim](https://github.com/brentp/hts-nim/) intended
to serve as examples for using `hts-nim` as well as to be useful tools.

These tools are:

```
hts-nim utility programs.
version: $version

	• bam-filter    : filter BAM/CRAM/SAM files with a simple expression language
	• count-reads   : count BAM/CRAM reads in regions given in a BED file
	• vcf-check     : check regions of a VCF against a background for missing chunks
```

each of these is described in more detail below.

# bam-filter

Use simple expressions to filter a BAM/CRAM file:

```
bam-filter

  Usage: bam-filter [options] <expression> <BAM-or-CRAM>

  -t --threads <threads>       number of BAM decompression threads [default: 0]
  -f --fasta <fasta>           fasta file for use with CRAM files [default: $env_fasta].
```

valid expressions may access the bam attibutes:

+  `mapq `/ `start `/ `pos `/ `end `/ `flag `/ `insert_size ` (where pos is the 1-based start)
+ `is_aligned` `is_read1` `is_read2` `is_supplementary` `is_secondary` `is_dup` `is_qcfail`
+ `is_reverse` `is_mate_reverse` `is_pair` `is_proper_pair` `is_mate_unmapped` `is_unmapped`

to use aux tags, indicate them prefixed with 'tag_', e.g.:

  tag_NM < 2. Any tag present in the bam can be used in this manner.

example:
```
bam-filter "tag_NM == 2 && tag_RG == 'SRR741410' && is_proper_pair" tests/HG02002.bam
```

# count-reads

Count reads reports the number of reads overlapping each interval in a BED file.

```
count-reads

  Usage: count-reads [options] <BED> <BAM-or-CRAM>

Arguments:                                                                                                                                                 

  <BED>          the bed file containing regions in which to count reads.
  <BAM-or-CRAM>  the alignment file for which to calculate depth.

Options:

  -t --threads <threads>      number of BAM decompression threads [default: 0]
  -f --fasta <fasta>          fasta file for use with CRAM files [default: ].
  -F --flag <FLAG>            exclude reads with any of the bits in FLAG set [default: 1796]
  -Q --mapq <mapq>            mapping quality threshold [default: 0]
  -h --help                     show help
```

This is output a line with a count of reads for each line in <BED>.

# vcf-check

`vcf-check` is useful as a quality control for large projects which have done variant calling in regions
where each region is called in parallel. With many regions, and large projects, some regions can error and
this might be unknown to the analyst.

This tools takes a background VCF, such as gnomad, that has full genome (though in some cases, users will
instead want whole exome) coverage and uses that as an expectation of variants. **If the background has many
variants across a long stretch of genome where the query VCF has no variation, we can expect that region is
missed in the query VCF.**

```
Check a VCF against a background to make sure that there are no large missing chunks.

  vcf-check

  Usage: vcf-check [options] <BACKGROUND_VCF> <VCF>

Arguments:                                                                                                                                                 
  <BACKGROUND_VCF>        population VCF/BCF with expected sites
  <VCF>                   query VCF/BCF to check

Options:

  -c --chunk <INT>        chunk size for genome [default: 100000]
  -m --maf <FLOAT>        allele frequency  cutoff [default: 0.1]
```

This will output a tab-delimited file of `chrom\tposition\tbackground-count\tquery-count`.

The user can find regions that might be problematic by plotting or with some simple `awk` commands.
