extern crate bamreadfilt;
extern crate byteorder;
extern crate rust_htslib;

use std::env;
use std::str;
use std::process;
use bamreadfilt::Config;
use bamreadfilt::BAMVCFRecord;
//use std::thread;
//use std::time::Duration;
//use byteorder::BigEndian;
//use byteorder::ReadBytesExt;


/*
   IMPORTANT!!!!

 * Another thing to consider:
 * We have pileup capabilities within rust (ie given a set of reads at an alignment pileup variants  or some shit)
 *   we could look into using this to skip mpileup directly or to rewrite it using our tool.
 *   but that sounds hard -- consider it! benchmark it!

 * the next logical thing to do is see if we can take a bam, return it as an iterator (if it isnt already), and then
 *  write that bam to disk.
 * 
 * the EXTRA information we will ultimately need is: base/ref at position of variation, qual at position of variation, cigar at position of variation, mate pair?
 */


fn main() {
    let args: Vec<String> = env::args().collect();

    // we will need to configure our own setup here.
    let config = Config::new(&args).unwrap_or_else(|err| {
             println!("Problem parsing arguments: {}", err);
             process::exit(1);
    });
    parse_bam(&config.bam, &config.vcf);

    // cool so we can successfully parse the AD field now. what do we have left to do with this part?
    //VCF gives us locations to go look for sites in the bam. next up is actually parsing the bam.
}

fn parse_bam(bamname: &str, vcf: &str){ 
    /*
     * Accept a bamname, set of vcf sites and their alt-refs (just in a form that we can filter for quickly, maybe even a struct)
     *    returns an iterator of parsed bam entries that are contained in vcfinfo
     *   
     *    requirements:
     *        should be able to apply filers onto said iterator as defined in our use cases
     *        should be easy to convert back into a bam protobuf or into a tfrecord protobuf
     */

    use rust_htslib::bam;
    use rust_htslib::bam::record::CigarString;
    use std::path::Path;
    let path = Path::new(bamname);
    let mut bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");

    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;
    let mut out = bam::Writer::from_path("test/out.bam", &head).expect("Error opening bam for writing.");


    let tid = 0;

    let beg: u32 = 10000;
    let end: u32 = 10075;


    // parse vcf -> vec
    use rust_htslib::bcf::Reader;
    let vcf_path = Path::new(vcf);
    let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); //originally said .ok().unwrap() which seems straight stupid
    let mut vcf_list: Vec<rust_htslib::bcf::Record> = vec![];
    for entry in vcf_itr.records() {
        // TODO turn this into a match and explicitly handle failure cases (by exiting?)
        vcf_list.push(entry.ok().unwrap());
    }

    // should probably figure out how to turn these into tests
    println!("hello? guys?");
    use rust_htslib::bam::record::Cigar;
    use rust_htslib::bam::Read;


    let mut bamvcfrecord = BAMVCFRecord::new(&mut bam, vcf_list);

    // we constructed bamvcfrecord to only allow things that contained variants!
    for (i, (item, vcf)) in bamvcfrecord.enumerate() {
        // all you need is read_pos (??? yay)
        if i % 1000 == 0 {
            println!("{} reads processed.", i);
        }
        let ref_pos = item.cigar().read_pos(vcf.pos, false, false);
        //ultimately we will still have to do soething fancy here.
        out.write(&item).unwrap();
    }


}


