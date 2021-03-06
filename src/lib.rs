extern crate byteorder;
extern crate rust_htslib;
extern crate argparse;

use argparse::{ArgumentParser, StoreTrue, Store, StoreFalse, StoreOption};
use std::collections::HashSet;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bcf;
use std::fs::File;


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
    pub stats: String,
    pub before_bed: String,
    pub after_bed: String,
    pub remove_duplicates: bool,
    pub raw_stats_fn: String,
    pub tag: String,
    pub max_sites: i64,
    pub only_alts: bool, 

}

impl Config {
    pub fn new() -> Config {
        let verbose = false;
        let bam_filename = "".to_string();
        let vcf_filename = "".to_string(); 
        let mapq = 0;
        let vbq = 0;
        let mrbq = 0.;
        let min_pir = 0;
        let max_pir = 10000000;
        let min_frag = 0;
        let max_frag = std::u32::MAX; // i would be this is getting filtered.
        let use_stdout: bool = false;
        let stats_before: bool = false;
        let stats_after: bool = false;
        let max_sites = -1;
        let bam_out_filename = "out.bam".to_string(); 
        let num_threads = 2;
        let stats  = "stats.txt".to_string();
        let before_bed = "before.bed".to_string();
        let after_bed  = "after.bed".to_string();
        let raw_stats_fn = "/dev/null".to_string();
        let tag = "".to_string();
        let remove_duplicates = false;
        let only_alts = false;
        let mut conf = Config { verbose, bam_filename, vcf_filename, mapq, vbq, mrbq, min_pir, max_pir, use_stdout, 
                                stats_before, stats_after, bam_out_filename, min_frag, max_frag, num_threads, stats, before_bed, after_bed, remove_duplicates, raw_stats_fn, tag, max_sites, only_alts };

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
            ap.refer(&mut conf.remove_duplicates)
                .add_option(&["--filter_duplicates", "--filtdup", "-r"], StoreTrue,
                "Filter all marked duplicates from the final bam. (Default: False)");
            ap.refer(&mut conf.only_alts)
                .add_option(&["--only_alts", "-A"], StoreTrue,
                "Filter any read which does not have the alt of interest. (Defaults: traffic-control: False, runway: True)");
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
            ap.refer(&mut conf.stats)
                .add_option(&["--stats"], Store,
                "File where we would like to place statistics.");
            ap.refer(&mut conf.before_bed)
                .add_option(&["--sites_before_bed"], Store,
                "Destination filename for sites present before filtering.");
            ap.refer(&mut conf.after_bed)
                .add_option(&["--sites_after_bed"], Store,
                "Destination filename for sites present after filtering.");
            ap.refer(&mut conf.raw_stats_fn)
                .add_option(&["--raw_stats_fn"], Store,
                "Destination filename for raw read-level statistics after filtering.");
            ap.refer(&mut conf.num_threads)
                .add_option(&["--threads", "-t"], Store,
                "Number of threads.");
            ap.refer(&mut conf.tag)
                .add_option(&["--keep-tag", "-k"], Store,
                "Include this tag as an output (requires raw_stats_fn)");
            ap.refer(&mut conf.max_sites)
                .add_option(&["--sample-k-sites"], Store,
                "Randomly downsamples vcfs to k entries");
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
            if alleles[0].len() > 1 {
                panic!("error, deletion detected in vcf. Insertions and deletions are not supported. {} {}", std::str::from_utf8(record.header().rid2name(chrm)).unwrap(), pos);
            } 
            if alleles[1].len() > 1 {
                panic!("error, insertion detected in vcf. Insertions and deletions are not supported. {} {}", std::str::from_utf8(record.header().rid2name(chrm)).unwrap(), pos);
            } 
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
    read_tracker: UniqueReadTree,
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
        //interesting, so what causes this to an hero?
        // this was inline with the end of the fetch call. .expect("Error fetching entry in bam -- check your bam for corruption.");     
        let fetched = match bam_reader.fetch(vcf_head.chrm, vcf_head.pos, vcf_head.pos + 1){
            Ok(val) => val,
            Err(e) => {eprintln!("CAUGHT AN ERROR IN NEW {}", e); panic!("error fetching");},
        };
        BAMVCFRecord { _bam: bam_reader, _vcf: vcf_reader, vcf_head: Some(vcf_head), read_tracker: UniqueReadTree::new(200) }

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

                //this was before debug route
                //self._bam.fetch(rid, entry.pos(), entry.pos() + 1).expect("Error fetching entry in bam -- check your bam for corruption.");

                //this is debug route
                let fetched = match self._bam.fetch(rid, entry.pos(), entry.pos() + 1){
                    Ok(val) => val,
                    Err(e) => {eprintln!("caught fetching error when advancing VCF {} {} {} {}", e, rid, entry.pos(), entry.pos()); panic!("error fetching");},
                };
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
    fn next(&mut self) -> Option<(bam::record::Record, VCFRecord, usize)>{ // We should make this a struct or wrap it in a Box type to get a constant size.
        // The other approach to this method would be to fetch all the reads at all sites, then remove non-unique reads (or to store the reads and remove duplicates as we go)
        //   Since we work with a sorted bam, the UniqueRead should be correct.
        let include_softclips = false;
        let include_dels = false;
        loop {
            match self.vcf_head {
                None => return None,
                Some(vcf_entry) => {
                    let next_record = self._bam.records().next();
                    match next_record {
                        None => { self._advance_vcf_head(); },
                        Some(val) => { 
                            match val { 
                                Ok(val) => {
                                    //  INSIDE our iterator. how can we move it such that it gets applied post-filter?
                                    // we are denoting position in the read as pir
                                    // some option for including variants that would occur in softmasked regions.
                                    let pir = match val.cigar().read_pos(vcf_entry.pos, include_softclips, include_dels).expect("Error parsing cigar string"){
                                        Some(val) => val as usize,
                                        None => { continue; }
                                    }; 
                                   
                                    let ur = UniqueRead::new(&val);
                                    if self.read_tracker.contains(&ur) {
                                        continue;

                                    // here is where I would make the change... basically we can either add a flag and use it here
                                    //      or we can pull this out and apply it as a filter. Filter seems smarter, but I believe
                                    
                                    } else{
                                        self.read_tracker.insert(UniqueReadItem::UniqueRead(ur));
                                        return Some((val, vcf_entry.clone(), pir)); 
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

/// mrbq returns the mean base quality of the passed `record`, beautiful.
pub fn mrbq(record: &bam::Record) -> f32 {
    let res: u32 = record.qual().into_iter().map( |&x| x as u32 ).sum(); 
    //return (res as f32) / (record.cigar_len() as f32)
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

#[derive(Clone,Eq,Hash,Debug,Ord,PartialEq,PartialOrd)]
struct UniqueRead {
    /// derived ords goes top to bottom, so by having pos first we basically say position is the most important thing for ordering, followed by qname and template position
    ///     if a read is the same as another, this is irrelevant
    ///     if a read has the same qname but has a different template position its teh read pair (and should have a different pos too)
    ///     if neither of those are true, then the position is the most important thing to consider for ordering as we will be using this to keep a size restricted BTree.
    chrm: i32,
    pos: i32,
    qname: String,
    template_pos: bool,
}

impl UniqueRead {
    fn new(item: &bam::Record) -> UniqueRead {
        let qname: String = String::from_utf8((&item.qname()).to_vec()).unwrap();
        UniqueRead{ chrm: item.tid(), pos: item.pos(), qname: qname, template_pos: item.is_first_in_template() }
    }
}

/// use this type for switching between bams and unique reads (useful for mocking bams for instance)
enum UniqueReadItem {
    BamRecord(bam::Record),
    UniqueRead(UniqueRead),
}

struct UniqueReadTree {
    tree: std::collections::HashSet<UniqueRead>,
    max_distance: i32,
}

impl UniqueReadTree {

    //  I haven't come up with a failure case so
    fn insert(&mut self, item: UniqueReadItem) {
        let unique_read = match item {
            UniqueReadItem::BamRecord(item) => {
                let qname: String = String::from_utf8((&item.qname()).to_vec()).unwrap();
                UniqueRead{ chrm: item.tid(), pos: item.pos(), qname: qname, template_pos: item.is_first_in_template() }
            },
            UniqueReadItem::UniqueRead(i) => i 
        };

        // there is certainly a faster way to do this, there should be an inner node that is not the smallest and < unique_read
        //      there is another solution that involves updating once we hit a capacity doing the exact same thing! then we dont update
        //      on every insertion, just on insertions where we have reached capacity although we would need to drop all elements in the tree < k somehow.
        loop {
            let smallest = match self.tree.iter().min() {
                Some(val) => 
                     val.clone(),
                None =>  unique_read.clone(),
            };
            //this should be a loop!
            if (unique_read.chrm == smallest.chrm) && ((unique_read.pos - smallest.pos ) <= self.max_distance ) {
                break;
            }
            else {
                self.tree.remove(&smallest); 
            }
        }

        //now that the tree has been pruned we add our read.
        self.tree.insert(unique_read);
    }

    /// constructor. choose max_distance to a size that is equal to readlength or the upper bound on read length.
    fn new(max_distance: i32) -> UniqueReadTree {
        return UniqueReadTree { tree: std::collections::HashSet::new(), max_distance: max_distance }
    } 

    fn iter(&self) -> std::collections::hash_set::Iter<UniqueRead> {
        return self.tree.iter()
    }

    fn len(&self) -> usize {
        return self.tree.len()
    }

    fn max(&self) -> Option<&UniqueRead> {
        return self.tree.iter().max()
    }

    fn min(&self) -> Option<&UniqueRead> {
        return self.tree.iter().min()
    }
    fn contains(&self, other: &UniqueRead ) -> bool {
        return self.tree.contains(other)
    }

}


//now we need a limited BTree type/function to handle this shit.
//basically we want a BTree with a max capacity where we effectively drop shit that is < K 
//      need to check the complexity of keepign this balanced
//      need to check the max memory usage of this kind of method


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

    #[test]
    fn unique_read_ord(){
        let ur_first = UniqueRead { chrm: 1, pos: 151, qname: "ayylmao".to_string(), template_pos: true};
        let ur_first_clone = UniqueRead { chrm: 1, pos: 151, qname: "ayylmao".to_string(), template_pos: true};
        let ur_second = UniqueRead { chrm: 1, pos: 152, qname: "ayylmao_diff".to_string(), template_pos: true };
        assert_ne!(ur_first, ur_second);
        assert_eq!(ur_first, ur_first_clone);
        assert_eq!(ur_first.cmp(&ur_second), std::cmp::Ordering::Less);
        assert_eq!(ur_second.cmp(&ur_first), std::cmp::Ordering::Greater);
    }

    #[test]
    fn unique_read_bst_update(){
        // not really a test. more just im working with it to figure out whats goig on.
        let mut tree:  UniqueReadTree  = UniqueReadTree::new(200);
        let ur_first = UniqueRead { chrm: 1, pos: 151, qname: "ayylmao".to_string(), template_pos: true};
        let ur_first_clone = UniqueRead { chrm: 1, pos: 151, qname: "ayylmao".to_string(), template_pos: true};
        let ur_first_pair = UniqueRead { chrm: 1, pos: 302, qname: "ayylmao".to_string(), template_pos: false};
        let ur_second = UniqueRead { chrm: 1, pos: 1502, qname: "ayylmao_diff".to_string(), template_pos: true };
        let ur_third = UniqueRead { chrm: 1, pos: 42000, qname: "ayylmao_diff".to_string(), template_pos: true };
        let ur_fourth = UniqueRead { chrm: 1, pos: 42069, qname: "ayylmao_diff".to_string(), template_pos: true };
        let ur_fifth = UniqueRead { chrm: 2, pos: 1, qname: "ayylmao_diffchrm".to_string(), template_pos: true };
        
        tree.insert(UniqueReadItem::UniqueRead(ur_first.clone())); 
        assert_eq!(tree.len(), 1);
        assert_eq!(tree.contains(&ur_first_clone), true);
        tree.insert(UniqueReadItem::UniqueRead(ur_first_clone.clone()));
        assert_eq!(tree.len(), 1);
        tree.insert(UniqueReadItem::UniqueRead(ur_first_pair.clone()));
        assert_eq!(tree.len(), 2);
        tree.insert(UniqueReadItem::UniqueRead(ur_third.clone()));
        assert_eq!(tree.len(), 1);
        assert_eq!(tree.contains(&ur_first_clone), false);
        tree.insert(UniqueReadItem::UniqueRead(ur_fourth.clone()));
        assert_eq!(tree.len(), 2);
        tree.insert(UniqueReadItem::UniqueRead(ur_fifth.clone())); // should drop all the others
        assert_eq!(tree.len(), 1);
        assert_eq!(tree.contains(&ur_third), false);
        assert_eq!(tree.contains(&ur_fourth), false);

        // this is nice with nocap
        for val in tree.iter(){
            println!("{:?}", val);
        }

    }

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
    fn running_average_update_range() {

        let mut ra = RunningAverage::new();
        for i in 0..1000000 {
            ra.update(i);
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
        let mut bamvcfrecord = BAMVCFRecord::new(&mut "test/test.bam", vcf_list);
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


