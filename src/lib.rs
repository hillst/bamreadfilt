#![feature(test)]
#![feature(rand)]

extern crate test;
extern crate rand;
extern crate byteorder;
extern crate rust_htslib;
extern crate argparse;

use argparse::{ArgumentParser, StoreTrue, Store, StoreFalse};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bcf;

#[derive(Debug,Clone)]
pub struct Config {
    pub verbose: bool,
    pub bam_filename: String , 
    pub vcf_filename: String , 
    pub mapq: u8,
    pub vbq: u8, 
    pub mrbq: f32, 
    pub min_pir: usize,
    pub max_pir: usize,
    pub use_stdout: bool,
    pub stats_before: bool,
    pub stats_after: bool, 
    pub bam_out_filename: String,
    pub min_frag: u32,
    pub max_frag: u32,
    pub num_threads: u8,

}

impl Config {
    pub fn new() -> Config {

        let mut verbose = false;
        let mut bam_filename = "".to_string();
        let mut vcf_filename = "".to_string(); 
        let mut mapq = 0;
        let mut vbq = 0;
        let mut mrbq = 0.;
        let mut min_pir = 0;
        let mut max_pir = 10000000;
        let mut min_frag = 0;
        let mut max_frag = 10000000;
        let mut use_stdout: bool = false;
        let mut stats_before: bool = true;
        let mut stats_after: bool = true;
        let mut bam_out_filename = "out.bam".to_string(); 
        let mut num_threads = 1;
        let mut conf = Config { verbose, bam_filename, vcf_filename, mapq, vbq, mrbq, min_pir, max_pir, use_stdout, 
                                stats_before, stats_after, bam_out_filename, min_frag, max_frag, num_threads };

        {  // this block limits scope of borrows by ap.refer() method
            let mut ap = ArgumentParser::new();
            ap.set_description("Welcome to bamreadfilt. This is a tool for filtering bams which contain variants. In the most simple use case, bamreadfilt simply returns
                               a new bam of reads that only contain known variants from a provided vcf. In addition, bamreadfilt provides a simple interface for gathering
                               statistics and filtering read by a variety of read-specific features.");
            ap.refer(&mut conf.verbose)
                .add_option(&["-v", "--verbose"], StoreTrue,
                "Be verbose");
            ap.refer(&mut conf.bam_filename)
                .add_option(&["--bam", "-b"], Store,
                "Filename of indexed bam for reading.")
                .required();
            ap.refer(&mut conf.vcf_filename)
                .add_option(&["--vcf"], Store,
                "Filename of vcf for reading.")
                .required();
            ap.refer(&mut conf.use_stdout)
                .add_option(&["--stdout", "-s", "-1"], StoreTrue,
                "Write bam to stdout insetead of to a file (to sort, for instance).");
            ap.refer(&mut conf.stats_before)
                .add_option(&["--stats_before"], StoreTrue,
                "Collect filter statistics before filtering.");
            ap.refer(&mut conf.stats_after)
                .add_option(&["--stats_after"], StoreTrue,
                "Collect filter statistics for reads after filtering.");
            //filters
            ap.refer(&mut conf.mapq)
                .add_option(&["--mapq"], Store,
                "Minimum MAPQ to keep read");
            ap.refer(&mut conf.vbq)
                .add_option(&["--vbq"], Store,
                "Minimum base-quality (Phred 33) score at the position of the variant to include a read.");
            ap.refer(&mut conf.mrbq)
                .add_option(&["--mrbq"], Store,
                "Minimum read base quality to include a read. Computes the mean of the full read, not just the mapped portion of the read.");
            ap.refer(&mut conf.min_pir)
                .add_option(&["--min_pir"], Store,
                "Minimum position of the variant in the read to include the read.");
            ap.refer(&mut conf.max_pir)
                .add_option(&["--max_pir"], Store,
                "Maximum fragment size of the read if paired.");
            ap.refer(&mut conf.min_frag)
                .add_option(&["--min_frag"], Store,
                "Minimum fragment size of the read if paired.");
            ap.refer(&mut conf.max_frag)
                .add_option(&["--max_frag"], Store,
                "Maximum position of the variant in the read to include the read.");
            ap.refer(&mut conf.num_threads)
                .add_option(&["--threads", "-t"], Store,
                "Number of threads.");
            ap.parse_args_or_exit();
        }
        conf

    }
}





/// Simple structure for capturing a current VCF entry (removes the INFO field effectively)
#[derive(Copy,Clone,PartialEq,Eq,Hash,Debug)]
pub struct VCFRecord {
    pub chrm: u32,
    pub pos: u32,
    pub ref_b: u8,
    pub alt_b: u8, 
}

impl VCFRecord {
    ///  Constructs a new record, panics if the reference ID is invalid
    pub fn new(record: &bcf::Record) -> VCFRecord {
        let chrm = match record.rid(){
            Some(val) => val,
            None => panic!("Error, missing reference ID"),
        };
        let pos = record.pos();
        let alleles = record.alleles(); // Vec<&[u8]> signature
        if alleles.len() > 2 {
            // this shold produce a recoverable error but uh I guess we wont worry about that for now
            panic!("error, detected non-diploid site");
        } else {
            let ref_b = alleles[0][0].clone();
            let alt_b = alleles[1][0].clone();
            VCFRecord { chrm: chrm, pos : pos, ref_b: ref_b, alt_b: alt_b } 
        }

    }

}

/// _bam: current location of the bam file descriptor
/// _vcf: Underlying vector of VCF records
/// vcf_head: Contains the current vcf entry to fetch records from
pub struct BAMVCFRecord {
    _bam: bam::IndexedReader,
    _vcf: Vec<bcf::Record>, 
    vcf_head: Option< VCFRecord >, //lol why is this option again? because in the last instance we want to trash it?
}


impl BAMVCFRecord {
    // two iteration rules
    //     no reads should be handled by the iterator, constructor should simply setup the start of our vcf
    // one: if there are no reads, get the next vcf entry and fetch reads (and check if they exist)
    // two: if there is a vcf entry and reads return the next read 

    pub fn new(bam_path: &str , mut vcf_reader: Vec<bcf::Record>) -> BAMVCFRecord {
        let vcf_head = match vcf_reader.pop(){
            Some(entry) => VCFRecord::new(&entry),
            None => panic!("No entries in VCF vector."), // is this really necessary? it is nice to exit first..
        };
        let mut bam_reader = bam::IndexedReader::from_path(bam_path).ok().expect("Error opening bam.");
        bam_reader.fetch(vcf_head.chrm, vcf_head.pos, vcf_head.pos + 1).expect("Error fetching entry in bam -- check your bam for corruption.");     //make sure this is right
        BAMVCFRecord { _bam: bam_reader, _vcf: vcf_reader, vcf_head: Some(vcf_head) }

    }

    /// updates the structs vcf_head to be the next entry. Also updates the _bam file descriptor.
    ///     We use this when the previous vcf_head has been exhausted (no more reads)
    ///     If we have consumed our full iterator, there are no more entries and we should return None to our iterator. 
    fn _advance_vcf_head(&mut self) {
        let new_head = match self._vcf.pop(){
            Some(entry) => {
                let rid = match entry.rid(){
                    Some(value) => value,
                    None => panic!("Invalid reference id in VCF"),
                };
                //update bam pointer
                self._bam.fetch(rid, entry.pos(), entry.pos() + 1).expect("Error fetching entry in bam -- check your bam for corruption.");
                Some(VCFRecord::new(&entry))
            },
            None => {
                None  // if the top of our vcf thing is ever head we know to exit!
            },
        };
        self.vcf_head = new_head;
    }



}


impl Iterator for BAMVCFRecord {
    type Item = (bam::record::Record, VCFRecord, usize);

    /// returns the current read (as a bam::record::Record), variant (as a VCFRecord), and position of the mutation in the read (usize).
    fn next(&mut self) -> Option<(bam::record::Record, VCFRecord, usize)>{ //eventually we are going to make this a fancier struct
        let include_softclips = false;
        let include_dels = false;
        loop {
            match self.vcf_head {
                None => return None,
                Some(vcf_entry) => {
                    match self._bam.records().next() {
                        None => { self._advance_vcf_head(); },
                        Some(val) => { 
                            match val { 
                                    Ok(val) => {
                                        // we are denoting position in the read as pir
                                        // some option for including variants that would occur in softmasked regions.
                                        let pir = match val.cigar().read_pos(vcf_entry.pos, include_softclips, include_dels).expect("Error parsing cigar string"){
                                            Some(val) => val as usize,
                                            None => { continue; }
                                        }; 

                                        // lol this is our default filter, i wonder if we can apply this somewhere else
                                        //      to be more idiomatic. for instance, if this our default filter, shouldnt it show up
                                        //      with the other filters rather tahn inside the iterator?
                                        if val.seq()[pir] == vcf_entry.alt_b {
                                            return Some((val, vcf_entry.clone(), pir)); 
                                        } else{
                                            continue;
                                        }
                                    }, //it may be wise to repack these with all relevant information in addition to bam
                                    Err(_) => panic!("Error reading bam, godspeed."),
                            }
                        }
                    }
                }
            }
        }
    }
}

pub fn lofreq_shitty() {
    /* How shitty is lofreq?
     *
     * (MRD)
     * Of all sites called by lofreq in the tumor (CA samples and cfDNA),
     *  how many (stats) were then found in the plasma? Of how many found in the plasma are noise?
     *
     *
     *
     *  we define noise as anything that fails a filter
     *
     *
     *

     list_of_ca_samples = [];
     for bam, vcf in ca_samples
         let bamvcf = bamvcfrecord(bam, vcf);

         // update_stats(x)
         //     i would like to collect the number of site found and the number of reads associated with each site

         bamvcf.map( |x| { stat_holder_before.update_stats(x) }
                    )
               .filter( |x| myfilters(x))
               .filter( .. )
               .filter( .. )
               .map(stat_holder_after.update_stats(x))

               bamreadfilt filtstats bams vcfs
    */ 


}

/// mrbq returns the mean base quality of the passed `record`
pub fn mrbq(record: &bam::Record) -> f32 {
    let res: u32 = record.qual().into_iter().map( |&x| x as u32 ).sum(); 
    return (res as f32) / (record.qual().len() as f32)
}


/// site_stats computes all stats of interest for a given dataset.
///     each metrics has a mean/variance
///     the counts are raw counts (hopefully we dont run out of memory!)
/// there is also a use-case where we may want site-specific statistics. But we have not addressed that :(
#[derive(Debug,Clone,Copy)]
pub struct SiteStats{ 
    //worthwhile to add something like median?
    //softmasking, mismatches, insertions?

    // getting number of sites not super intuitive, probably should be its own stat?
    //      just a collection of sites? does it need its own type? we really just want a 'set-like' object. 

    pub num_sites: u64, 
    pub num_reads: u64, 
    pub mrbq: RunningAverage,
    pub vbq: RunningAverage,
    pub mapq: RunningAverage,
    pub pir: RunningAverage,   
    pub duplicated: RunningAverage,
    pub paired: RunningAverage,
    pub frag_size: RunningAverage,
}

//usage is still map, filter, filter, map
//so ultimately we are 'mapping' chained update calls

impl SiteStats{

    pub fn new() -> SiteStats {
        return SiteStats { 
                    num_sites: 0, num_reads: 0, mrbq: RunningAverage::new(), vbq: RunningAverage::new(),
                    mapq: RunningAverage::new(), pir: RunningAverage::new(), duplicated: RunningAverage::new(),
                    paired: RunningAverage::new(), frag_size: RunningAverage::new(),
               }
    }

    //for (i, (item, _vcf, _pir)) in bamvcfrecord
    pub fn collect_stats_stream(&mut self, item: &bam::Record, vcf: &VCFRecord, pir: usize) {
        self.num_reads += 1; //if unique..... check to make sure this is the same as the number of reads we get back.
        self.mrbq.update(mrbq(item) as f64);
        self.vbq.update( item.qual()[pir] as f64);
        self.mapq.update( item.mapq() as f64);
        self.pir.update( pir as f64 );
        self.duplicated.update( item.is_duplicate() as u32 as f64 );
        self.paired.update( item.is_paired() as u32 as f64 );
        if item.is_paired() {
            self.frag_size.update( item.insert_size().abs() as f64 );
        }
    }
}

#[derive(Debug,Clone,Copy)]
pub struct RunningAverage {
    count: u64,
    mean: f64,
    m2: f64,
}

/// computes the moments online of a stream. 
/// # example
/// 
/// let mut stats = SiteStats::new();
/// stuff.map( | read, vcf, pir | { stats.mapq.update(read.mapq)} 
///     ...
///     .collect()
/// let (mapq_mean, mapq_var) = stats.mapq.finalize();
/// 
impl RunningAverage {
    pub fn new() -> RunningAverage {
        return RunningAverage{ count:0, mean:0.0, m2:0.0 }        
    }

    /// Updates the running average -- interestingly, we might want a generic type for this (int vs floats?)
    /// algorithm taken from wikipedia page, the numerically stable version was selection to remove subtraction of very small numbers (precision issue)
    pub fn update(&mut self, new_value: f64){
        self.count += 1; //this guarantees count is never 0 correct?
        let delta = new_value as f64 - self.mean;              // -1
        self.mean = self.mean + (delta / self.count as f64); 
        let delta2 = new_value as f64 - self.mean;             
        self.m2 = self.m2 + delta * delta2;                    

    }

    /// returns the mean and variance for this running average. if the number of samples is less than 2 (ie no variance) then we return None.
    ///     use this function when you are done updating the struct :)
    pub fn finalize(&self) -> Option<(f64, f64)>{
        if self.count < 2 {
            None
        } else {
            let variance = self.m2 / (self.count) as f64; // (self.count - 1 was here before...)
            return Some((self.mean, variance))
        }
    }

}

pub fn tuple_to_csv(tuple: (f64, f64)) -> String{
    format!("{},{}", tuple.0, tuple.1)

}


#[cfg(test)]
mod tests {
    //need to come up with a test suite for our main functionality (that or its all encapsulated in the library we are using)
    //      we dont actually do much that could be 'buggy' since there isnt math or anything occuring.
    //      could/should just make sure we correctly get 'None' and panics when we should ?
    use super::*;
    use test::Bencher;

    /// tests the range ( 1, 6 )
    #[test]
    fn running_average_update_simple() {
        let mut ra = RunningAverage::new();
        ra.update(2.);
        assert_eq!(ra.mean, 2.0);
        ra.update(1.); 
        assert_eq!(ra.mean, 1.5);
        ra.update(3.);
        assert_eq!(ra.mean, 2.0);
        ra.update(4.); 
        assert_eq!(ra.mean, 2.5);
        ra.update(5.); 
        let (mean, var) = match ra.finalize(){
            Some((m, v)) => (m, v),
            None => {panic!("error, got None type from finalize during testing");}
        };
        println!("mean: {}, var: {}", mean, var);
        assert_eq!(mean, 3.0);
        assert_eq!(var, 2.0);
    }

    #[test]
    fn running_average_update_rand() {
        use rand::Rng;

        let mut ra = RunningAverage::new();
        let mut rng = rand::thread_rng();
        for i in 0..1000000 {
            ra.update(rng.gen::<u32>() as f64);
        }
        let (mean, var) = match ra.finalize(){
            Some((m, v)) => (m, v),
            None => {panic!("error, got None type from finalize during testing");}
        };

        assert_ne!(mean, 0.);
        assert_ne!(var, 0.);
    }

    #[test]
    #[should_panic]
    fn running_average_update_empty() {
        let mut ra = RunningAverage::new();
        let (mean, var) = match ra.finalize(){
            Some((m, v)) => (m, v),
            None => {panic!("error, got None type from finalize during testing");}
        };
    }

    #[test]
    fn running_average_update_constant() {
        let mut ra = RunningAverage::new();
        ra.update(2.);
        assert_eq!(ra.mean, 2.0);
        ra.update(2.); 
        assert_eq!(ra.mean, 2.0);
        ra.update(2.);
        assert_eq!(ra.mean, 2.0);
        ra.update(2.); 
        assert_eq!(ra.mean, 2.0);
        ra.update(2.); 
        let (mean, var) = match ra.finalize(){
            Some((m, v)) => (m, v),
            None => {panic!("error, got None type from finalize during testing");}
        };
        println!("mean: {}, var: {}", mean, var);
        assert_eq!(mean, 2.0);
        assert_eq!(var, 0.0);

    }


    #[test]
    fn stress_running_average_update() {
        let mut ra = RunningAverage::new();
        for i in 0..1000000 { // 1 BILLION iterations
            ra.update(2.);
        }
        let (mean, var) = match ra.finalize(){
            Some((m, v)) => (m, v),
            None => {panic!("error, got None type from finalize during testing");}
        };
        println!("mean: {}, var: {}", mean, var);
        assert_eq!(mean, 2.0);
        assert_eq!(var, 0.0);
    }

    //future benchmarks should include things like:
    //      vcfs with different number of sites (basically dense)
    //      how frequent are reads containing mutations (reads with muts / reads)
    //      how dense are reads containing mutations -- coverage (reads with muts / muts) 
    //      different number of reads (does increasing bamsize but fixing other variables keep performance?)
    #[bench]
    fn iterate_bamvcf(b: &mut Bencher) {
        // we are supposed to provide a closure to b.iter() that does one step of 'work'

        println!("Starting benchmark test:");
        use rust_htslib::bam;
        use std::path::Path;

        //prep bam
        let path = Path::new("test/test.bam");
        let mut bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");

        //prep vcf
        use rust_htslib::bcf::Reader;
        let vcf_path = Path::new("test/test4k.vcf");
        let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); //originally said .ok().unwrap() which seems straight stupid
        let mut vcf_list: Vec<rust_htslib::bcf::Record> = vec![];
        for entry in vcf_itr.records() {
            vcf_list.push(entry.ok().unwrap());
        }

        //prep iterator
        let mut bamvcfrecord = BAMVCFRecord::new(&mut bam, vcf_list);
        let mut test = 0;
        b.iter( || {
            match bamvcfrecord.next(){
                Some((ref item, vcf, pir)) => {
                    test += pir;
                },
                None => {return;},
            };
        });
    }
}


