# bamreadfilt

Bamreadfilt is a library/software suite for filtering BAM (https://samtools.github.io/) files according to read-specific filters related to variants. 

## Design goals

We are building this tool so there is a standard application and commandline interface for analyzing bams emphasizing specific variations. 

- Be idiomatic
- Be safe
- Support typical bioinformatics ecosystems (ie stdin <-> bams)
- Write strong tests
- Be parallel

## Building

```bash
git clone <thisrepo>
cd bamreadfilt/
cargo build --release
```
## Usage

The most basic usecase of bam read filt is to filter a bam for *only* reads that contain a variant of interest. Each use of bamreadfilt requires two things, a bam to be filtered, and a vcf/bed containing variants of interest. 

```bash
./bamreadfilt --bam test.bam --vcf test.vcf
```

That's it! Two files are written:
- stats.txt 
  This contains statistics (mean/variance) about the available filters and about the number of sites and reads before/after.
- out.bam
  This is a bam format corresponding to unfiltered reads.


## Filtering by thresholds

In the default case, brf (bamreadfilt) filters for reads that contain a target variation. However, we also have several additional methods to use when filtering our bam. These are:

```
  --mapq MAPQ           Minimum MAPQ to keep read
  --vbq VBQ             Minimum base-quality (Phred 33) score at the position
                        of the variant to include a read.
  --mrbq MRBQ           Minimum read base quality to include a read. Computes
                        the mean of the full read, not just the mapped portion
                        of the read.
  --min_pir MIN_PIR     Minimum position of the variant in the read to include
                        the read.
  --max_pir MAX_PIR     Maximum fragment size of the read if paired.
  --min_frag MIN_FRAG   Minimum fragment size of the read if paired.
  --max_frag MAX_FRAG   Maximum position of the variant in the read to include
                        the read.
 ```
Each one of these filters applies an additional filter to the resulting bam. 

For instance, if a user wants to filter by variant base quality >= 20 and MAPQ >= 40, you would use:

```bash
./bamreadfilt --bam test.bam --vcf test.vcf --mapq 40 --vbq 20
```

The result would be written to test/new_out.bam and stats.txt respectively. Unfortunately however, the resulting bam is not sorted (and not indexed!), to receive a sorted bam you can write to standard out with the --stdout flag. Unfortunately using this with samtools sort  **does not play work** with a limited amount of system memory. Avoid using --stdout to sort your result.

While the tool works great for bam/vcf pairs, there is no way currently to process more than one bam at a time. This could lead to frustration when collecting statistics. Fortunately there is a flag for a custom statistics file!

```bash
./bamreadfilt --bam test.bam --vcf test.vcf --mapq 40 --vbq 20  --stats mystats.txt 
```

Now results are written to mystats.txt instead of to the default stats.txt.

## Multithreading

To use bamreadfilt with threads, simply add the -t argument specifying a number of threads. The filters are not distributed and we are *guaranteed* there is no race condition when computing the statistics.

```bash 
./bamreadfilt --bam test.bam --vcf test.vcf --mapq 40 --vbq 20  -t 8 
```
