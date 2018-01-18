extern crate bamreadfilt;
extern crate byteorder;
extern crate rust_htslib;

use std::collections::HashSet;
use bamreadfilt::Config;
use bamreadfilt::tuple_to_csv;
use bamreadfilt::BAMVCFRecord;

fn main() {
    let config = Config::new();
    run(config);
}

fn run(config: Config){ 
    /*
     * Accept a bamname, set of vcf sites and their alt-refs (just in a form that we can filter for quickly, maybe even a struct)
     *    returns an iterator of parsed bam entries that are contained in vcfinfo
     *   
     *    requirements:
     *        should be able to apply filers onto said iterator as defined in our use cases
     *        should be easy to convert back into a bam protobuf or into a tfrecord protobuf
     */

    use rust_htslib::bam;
    use std::path::Path;

    use rust_htslib::bcf::Reader;
    use rust_htslib::bam::Read;

    let path = Path::new(&config.bam_filename);
    let mut bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");

    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;

    //need something here to determine between stdout (i think we can get a file handle for that) and an actual path.
    //  ideally our config struct will handle this.

    let mut out = match config.use_stdout {
        false => bam::Writer::from_path("test/new_out.bam", &head).expect("Error opening bam for writing."), 
        true => bam::Writer::from_stdout(&head).expect("Error opening stdout for writing."),
    };

    // parse vcf -> vec
    let vcf_path = Path::new(&config.vcf_filename);
    let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); 
    let mut vcf_list: Vec<rust_htslib::bcf::Record> = vec![];
    for entry in vcf_itr.records() {
        // TODO turn this into a match and explicitly handle failure cases (by exiting?)
        //      I guess it just isnt clear to me how this fails.
        vcf_list.push(entry.ok().unwrap());
    }

    let bamvcfrecord = BAMVCFRecord::new(&mut bam, vcf_list);

    let mut before_stats = bamreadfilt::SiteStats::new();
    let mut after_stats = bamreadfilt::SiteStats::new();

    let mut before_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();
    let mut after_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();


    //DAMN this is fucking beautiful -- need to paramterize it though.
    for (i, (item, _vcf, _pir)) in bamvcfrecord
                                        .map(   | (_i, _v, p) | { 
                                            before_stats.collect_stats_stream(&_i, &_v, p);
                                            before_site_stats.insert(_v.clone());
                                            (_i, _v, p)
                                        })
                                        .filter(| &(ref _i, _v, p)| { p >= config.min_pir && p <= config.max_pir}) //pir filter
                                        .filter(| &(ref i, _v, _p)| { 
                                            if i.is_paired() {
                                                i.insert_size().abs() as u32  >= config.min_frag && i.insert_size().abs() as u32 <= config.max_frag

                                            } else {
                                                true
                                            }
                                        }) //pir filter
                                        .filter(| &(ref i, _v, p) | {                     //vbq filter
                                            i.qual()[p] >= config.vbq
                                        })
                                        .filter(| &(ref i, _v, _p)| { i.mapq() >= config.mapq })    //mapq filter
                                        .filter(| &(ref i, _v, _p)| { bamreadfilt::mrbq(i) >= config.mrbq }) //mrbq filter
                                        .map(   | (_i, _v, p) | { 
                                            after_stats.collect_stats_stream(&_i, &_v, p);
                                            after_site_stats.insert(_v.clone());
                                            (_i, _v, p)
                                        })
                                        .enumerate() 
    {
        if config.verbose && i % 1000 == 0 {
            eprintln!("{} reads processed.", i);
        }
        out.write(&item).unwrap(); // writes filtered reads.
    }


    write_statistics(&"stats.txt".to_string(), before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
}

//  probably should put the before/after sites somewhere proper.
//  would be nice if we could write the config somewhere too! Then we would have a list of applied filters in addition to before/after stats.
pub fn write_statistics(output_filename: &String, before_stats: bamreadfilt::SiteStats, after_stats: bamreadfilt::SiteStats, before_sites: u64, after_sites: u64) -> Result<(),std::io::Error> {
    use std::io::Write;
    use std::fs::File;


    let mut f = File::create(output_filename)?;
    writeln!(&mut f, "setting,mrbq_mean,mrbq_var,vbq_mean,vbq_var,mapq_mean,mapq_var,duplicated_mean,duplicated_var,frag_size,num_sites,num_reads").unwrap();
    writeln!(&mut f, "before,{},{},{},{},{},{},{}", 
               tuple_to_csv(before_stats.mrbq.finalize().unwrap()), 
               tuple_to_csv(before_stats.vbq.finalize().unwrap()), 
               tuple_to_csv(before_stats.mapq.finalize().unwrap()), 
               tuple_to_csv(before_stats.duplicated.finalize().unwrap()),
               tuple_to_csv(before_stats.frag_size.finalize().unwrap()),
               before_sites,
               before_stats.num_reads
           )?;
    writeln!(&mut f, "after,{},{},{},{},{},{},{}", 
               tuple_to_csv(after_stats.mrbq.finalize().unwrap()), 
               tuple_to_csv(after_stats.vbq.finalize().unwrap()), 
               tuple_to_csv(after_stats.mapq.finalize().unwrap()), 
               tuple_to_csv(after_stats.duplicated.finalize().unwrap()),
               tuple_to_csv(after_stats.frag_size.finalize().unwrap()),
               after_sites,
               after_stats.num_reads
           )?;
    Ok(())

}


