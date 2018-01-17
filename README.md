# bamreadfilt

Bamreadfilt is a library/software suite for filtering BAM (https://samtools.github.io/) files according to read-specific filters related to variants. 

## Design goals

We are building this tool so there is a standard application and commandline interface for analyzing bams emphasizing specific variations. 

- Be idiomatic
- Be safe
- Support typical bioinformatics ecosystems (ie stdin <-> bams)
- Write strong tests
- Strong design (use driven development!)


## Building

```bash
git clone <thisrepo>
cd bamreadfilt/
cargo build --release
```
## Usage

The most basic usecase of bam read filt is to filter a bam for *only* reads that contain a variant of interest. Each use of bamreadfilt requires two things, a bam to be filtered, and a vcf/bed containing variants of interest. 

```bash
target/release/bamreadfilt test.bam test.vcf
```

## Additional Use Cases
- enforce uniqueness of reads
- collect statistics on vcf entries (number of reads per site, how many sites were present)
- collect statistics on reads (variant base qualitity, mean read base quality, mapping quality... etc)
- collect statistics for each site in the vcf as above (?)
- compute statistics for a single bam and vcf
- compute statistics on multiple bams and vcfs

- A user has a bam and would like to only select reads that have a base quality score > 20 at the position of variants.
- Same as previous, but would like to filter for MAPQ < 40
- only include reads which have no softclipping
- apply the VBQ and MQ filter to the result of a mpileup call
- tfrecords ?

