# Runway and Traffic Control

Runway and Traffic control are a library/software suite for filtering BAM (https://samtools.github.io/) files according to read-specific filters related to variants. 
traffic-control was designed for re-analyzing genotypes by fetching all reads which map to particular variants.
runway was designed for analyzing mutations. This is done by fetching reads which harbor a specific alternate base (alt).

## Design goals

We are building this tool so there is a standard application and commandline interface for analyzing bams emphasizing specific variations. 

- Be idiomatic
- Be safe
- Support typical bioinformatics ecosystems (ie stdin <-> bams)
- Write strong tests
- Be parallel

## Building

Built with cargo 0.23.0-nightly (e447ac7e9 2017-09-27)

```bash
git clone <thisrepo>
cd bamreadfilt/
cargo build --release
target/release/runway --help
```
## Usage

The most basic usecase of bam read filt is to filter a bam for *only* reads that contain a variant of interest. Each use of bamreadfilt requires two things, a bam to be filtered, and a vcf/bed containing variants of interest. 

```bash
./runway -t 2 --bam test.bam --vcf test.vcf
```

That's it! Two files are written:
- stats.txt 
  This contains statistics (mean/variance) about the available filters and about the number of sites and reads before/after.
- out.bam
  This is a bam format corresponding to unfiltered reads.


## Filtering by thresholds

In the default case, runway (formerly bamreadfilt) filters for reads that contain a target variation. However, we also have several additional methods to use when filtering our bam. Use --help
to see a complete list of filtering criterias.

For instance, if a user wants to filter by variant base quality >= 20 and MAPQ >= 40, and only variants that occur in the first 150 bases you would use:

```bash
./runway --bam test.bam --vcf test.vcf --mapq 40 --vbq 20 --max_pir 150
```

The result would be written to test/new_out.bam and stats.txt respectively. Unfortunately however, the resulting bam is not sorted (and not indexed!), to receive a sorted bam you can write to standard out with the --stdout flag. Unfortunately using this with samtools sort  **does not play work** with a limited amount of system memory. Avoid using --stdout to sort your result.

Runway also allows users to collect statistics about each variant after filtering.
```bash
./runway --bam test.bam --vcf test.vcf --mapq 40 --vbq 20  --stats mystats.txt 
```

## Multithreading

Runway and traffic-control only will execute with threads, by default we use 2 threads, but more can be specified. The example below executes runway with 8 threads.

```bash 
./runway --bam test.bam --vcf test.vcf --mapq 40 --vbq 20  -t 8 
```
