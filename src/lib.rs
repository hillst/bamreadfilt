#![feature(test)]
extern crate test;
extern crate byteorder;
extern crate rust_htslib;

use std::error::Error;
use std::fs::File;
use std::io::prelude::*;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bcf;
use std::path::Path;




pub struct Config {
    pub bam: String,
    pub vcf:String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str>{
        if args.len() < 3 {
            return Err("not enough arguments"); // why did I need return + statement here instead of just expr
        }
        let bam = args[1].clone();
        let vcf = args[2].clone();
        Ok(Config { bam, vcf})
        
    }
}


#[derive(Copy,Clone)]
pub struct VCFRecord {
    pub chrm: u32,
    pub pos: u32,
    pub ref_b: u8,
    pub alt_b: u8, 
}

impl VCFRecord {
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


pub struct BAMVCFRecord<'a> {
    //something about lifetimes!!! how do we access our current 'thing'
    _bam: &'a mut bam::IndexedReader,
    _vcf: Vec<bcf::Record>, 
    vcf_head: Option< VCFRecord >, //lol why is this option again? because in the last instance we want to trash it?
    _EOL: bool,
}


impl <'a > BAMVCFRecord<'a> {
    // two iteration rules
    //     no reads should be handled by the iterator, constructor should simply setup the start of our vcf
    // one: if there are no reads, get the next vcf entry and fetch reads (and check if they exist)
    // two: if there is a vcf entry and reads return the next read 
    pub fn new(bam_reader: &'a mut bam::IndexedReader, mut vcf_reader: Vec<bcf::Record>) -> BAMVCFRecord <'a> { // if i use a reference will i need a lifetime?
        //get first vcf entry
        //set iterator head using fetch
        let vcf_head = match vcf_reader.pop(){
            Some(entry) => VCFRecord::new(&entry),
            None => panic!("No entries in VCF vector."), // is this really necessary? it is nice to exit first..
        };

        bam_reader.fetch(vcf_head.chrm, vcf_head.pos, vcf_head.pos + 1);     //make sure this is right
        BAMVCFRecord { _bam: bam_reader, _vcf: vcf_reader, vcf_head: Some(vcf_head), _EOL: false }
    }


    fn _advancevcf_head(&mut self) {
        let new_head = match self._vcf.pop(){
            Some(entry) => {
                let rid = match entry.rid(){
                    Some(value) => value,
                    None => panic!("Invalid reference id in VCF"),
                };
                //update bam pointer
                self._bam.fetch(rid, entry.pos(), entry.pos() + 1);
                Some(VCFRecord::new(&entry))
            },
            None => {
                None  // if the top of our vcf thing is ever head we know to exit!
            },
        };
        self.vcf_head = new_head;
    }



}


fn _contains_variant(bam_read: &bam::record::Record, vcf: &VCFRecord) -> Option<u32> {
    /// Returns the position in bam_read that contains vcf_entry. Otherwise returns None.

    /*
     * bam.pos, for entry (something?) in bam, walk matched bases until pos = vcf.pos, 
     *   return if alt = vcf.alt, otherwise None
     * vcf.pos
     *
     */ 


    None

}


impl <'a> Iterator for BAMVCFRecord<'a> {
    /// iterator that seeks through a bam in the order resolved by a vcf record


    /*
     * The next thing this piece of code needs to do is....
     *      either figure out how to dervice copy/clone traits for bcf::record
     *        or manually parse and create our own vcf struct (i personally like this option more since we want to use this stuff downstream)
     */
    type Item = (bam::record::Record, VCFRecord);

    fn next(&mut self) -> Option<(bam::record::Record, VCFRecord)>{ //eventually we are going to make this a fancier struct
        loop {
            match self.vcf_head {
                None => return None,
                Some(vcf_entry) => {
                    match self._bam.records().next() {
                        None => { self._advancevcf_head(); },
                        Some(val) => { 
                            match val { 
                                    Ok(val) => {
                                        // we are denoting position in the read as pir
                                        // some option for including variants that would occur in softmasked regions.
                                        let pir = match val.cigar().read_pos(vcf_entry.pos, false, false).expect("Error parsing cigar string"){
                                            Some(val) => val as usize,
                                            //None means this read didnt have the annotated variant!
                                            None => {continue;}
                                        }; 

                                        if val.seq()[pir] == vcf_entry.alt_b {
                                            return Some((val, vcf_entry.clone())); 
                                        } else{
                                            continue;
                                        }
                                    }, //it may be wise to repack these with all relevant information in addition to bam
                                    Err(val) => panic!("Error reading bam, godspeed."),
                            }
                        }
                    }
                }
            }
        }
    }
}


pub fn run() {
    println!("{}", "hello world");
}


#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;
    #[test]
    fn one_result() {

    }

    #[bench]
    fn iterate_bamvcf(b: &mut Bencher) {
        // we are supposed to provide a closure to b.iter() that does one step of 'work'

        println!("Starting benchmark test:");
        use rust_htslib::bam;
        use rust_htslib::bam::record::CigarString;
        use std::path::Path;

        //prep bam
        let path = Path::new("test/test.bam");
        let mut bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");

        //prep vcf
        use rust_htslib::bcf::Reader;
        let vcf_path = Path::new("test/test.vcf");
        let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); //originally said .ok().unwrap() which seems straight stupid
        let mut vcf_list: Vec<rust_htslib::bcf::Record> = vec![];
        for entry in vcf_itr.records() {
            vcf_list.push(entry.ok().unwrap());
        }

        //prep iterator
        use rust_htslib::bam::record::Cigar;
        use rust_htslib::bam::Read;
        let mut bamvcfrecord = BAMVCFRecord::new(&mut bam, vcf_list);

        b.iter( || {
            match bamvcfrecord.next(){
                Some((item, vcf)) => {
                    item.cigar().read_pos(vcf.pos, false, false)
                },
                None => {return;},
            };
        });


    }

}
