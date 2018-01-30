extern crate bamreadfilt;
extern crate byteorder;
extern crate rust_htslib;

use std::collections::HashSet;
use bamreadfilt::Config;
use bamreadfilt::tuple_to_csv;
use bamreadfilt::BAMVCFRecord;

fn main() {
    let config = Config::new();

    let num_threads = config.num_threads;
    if num_threads > 1 {
        mscp_run(config, num_threads as usize);
    } else{
        run(config);
    }
}

/// runs our config in a mscp setting. this also works for a single-threaded job (n_thread = 1).
///     currently the only work that is 'parallel' is identifying reads which contain a specific variant. we should ultimately look to move the filters/stats there as well.
///     this kind of job would be best if we had a server/receiver setup going on.
///         it also seems like currently it just waits for everything to finish before trying to consume -- it would be better if it didnt block :)
fn mscp_run(config: Config, n_threads: usize){
    use std::thread;
    use std::sync::mpsc;
    use std::path::Path;
    use rust_htslib::bcf::Reader;
    use rust_htslib::bam::Read;
    use rust_htslib::bam;
    use rust_htslib::bam::record::Aux;
    let bam_filename = config.bam_filename.clone();
    
    //TODO duplicated code
    // parse vcf -> vec
    let vcf_path = Path::new(&config.vcf_filename);
    let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); 

    let path = Path::new(&bam_filename);//why is it only using move semantics lol
    let bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");
    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;

    let mut out = match config.use_stdout {
        false => bam::Writer::from_path("new_out.bam", &head).expect("Error opening bam for writing."), 
        true => bam::Writer::from_stdout(&head).expect("Error opening stdout for writing."),
    };
    
    let mut before_stats = bamreadfilt::SiteStats::new();
    let mut after_stats = bamreadfilt::SiteStats::new();

    let mut before_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();
    let mut after_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();

    let mut vec_list: Vec<Vec<rust_htslib::bcf::Record>> = vec![];
    for _ in 0..n_threads{
        vec_list.push(vec![]);
    }


    for (i, entry) in vcf_itr.records().enumerate() {
        // I'm almost certain we will want these to be continuous in memory or atleast adjacent
        vec_list[i % n_threads].push(entry.ok().unwrap());
    }

    if vec_list.len() < n_threads {
        panic!("Number of VCF entries < number of threads. Is this VCF okay?");
    }

    let mut handles = Vec::with_capacity(n_threads);

    //THREADS start here.

    let (tx, rx) = mpsc::channel();
    for _ in 0..n_threads {
        let vcf_list = vec_list.pop().unwrap();
        let bam_filename = bam_filename.clone();
        let tx = tx.clone();
        handles.push(thread::spawn(move || {
            let bamvcfrecord = BAMVCFRecord::new(&bam_filename, vcf_list);
            for (item, _vcf, _pir) in bamvcfrecord
            // we can alternatively have filters apply here, but then we will have trouble with computingstatistics.
            {
                tx.send((item, _vcf, _pir)).unwrap();
            }
        }));
    }


    // this part feels like something that is easy to modify and something we could bake into bamvcf iterator
    let config = config.clone();
    let consumer = thread::spawn(move || {
        for (i, (mut item, _vcf, _pir)) in rx.iter()
            .map(   | (_i, _v, p) | { 
                before_stats.collect_stats_stream(&_i, &_v, p); //NOTE this works because it is the consumer, so we have no race conditions at this point
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
            .enumerate() {

            if config.verbose && i % 1000 == 0 {
                eprintln!("{} reads processed.", i);
            }
            item.push_aux( "B0".as_bytes(), &Aux::Integer(_vcf.pos as i64) );
            item.push_aux( "B1".as_bytes(), &Aux::Char(_vcf.alt_b) );
            item.push_aux( "B2".as_bytes(), &Aux::Integer(_pir as i64) );
            out.write(&item).unwrap(); // writes filtered reads.
        }
            
        eprintln!("Done.");
        //write_statistics(&"stats.txt".to_string(), before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
        write_statistics(&config.stats, before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
        write_bed(&config.before_bed, &before_site_stats).unwrap();
        write_bed(&config.after_bed, &after_site_stats).unwrap();

    });

    for handle in handles {
        handle.join().unwrap();
    }
    drop(tx); // do not remove this or you WILL DEADLOCK
    consumer.join().unwrap();

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
    use rust_htslib::bam::record::Aux;

    let bam_filename = config.bam_filename.clone();
    //only used for building header
    let path = Path::new(&bam_filename); //why is it only using move semantics lol
    let bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");
    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;

    let mut out = match config.use_stdout {
        false => bam::Writer::from_path("new_out.bam", &head).expect("Error opening bam for writing."), 
        true => bam::Writer::from_stdout(&head).expect("Error opening stdout for writing."),
    };

    // parse vcf -> vec
    let vcf_path = Path::new(&config.vcf_filename);
    let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); 
    let mut vcf_list: Vec<rust_htslib::bcf::Record> = vec![];
    for entry in vcf_itr.records() {
        // TODO turn this into a match and explicitly handle failure cases (by exiting?)
        vcf_list.push(entry.ok().unwrap());
    }

    let bamvcfrecord = BAMVCFRecord::new(&bam_filename, vcf_list);

    let mut before_stats = bamreadfilt::SiteStats::new();
    let mut after_stats = bamreadfilt::SiteStats::new();

    let mut before_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();
    let mut after_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();


    for (i, (mut item, _vcf, _pir)) in bamvcfrecord
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
            //use std::mem::size_of;
            
            //println!( "size of {}", size_of::<(bam::record::Record, bamreadfilt::VCFRecord, usize)>() );
        }
        // these work so we can let them just kinda exist
        //
        // println!( "{:?}", item.seq().as_bytes() );
        // println!( "{:?}", item.cigar());
        // println!( "{:?}", item.qual());
        //
        // This also works
        //   
        //
        item.push_aux( "B0".as_bytes(), &Aux::Integer(_vcf.pos as i64) );
        item.push_aux( "B1".as_bytes(), &Aux::Char(_vcf.alt_b) );
        item.push_aux( "B2".as_bytes(), &Aux::Integer(_pir as i64) );


        // so this approach would be along the lines of just going from the bam and using the read tags to get the pir *directly*
        //      if we do it this way, we dont havet o fumble computing pir ourselves but we still would need to break for the different cigar strings.

        // the alternate is to use message passing by serializing into json and then to deserialize in python and to just roll with it. 
        //      I still kind of like this more? i dont really know.. I'm going to sleep on it.

        //use bamreadfilt::tfrecord;

        // we would want to expand this out, or create a way to iterate between these two thigns, or align them. these are what need to be aligned.
        //      basically we just need to pay attention to deletions

        // after working this out more (since there is no reference sequence here) it seems like we are going to want to use the other approach (keep mine and use work 
        //              with tags)


        //println!("{:?} {:?}", item.cigar().len(), item.seq().len() );

        out.write(&item).unwrap(); // writes filtered reads.

    }

    write_statistics(&config.stats.to_string(), before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
    write_bed(&config.before_bed, &before_site_stats).unwrap();
    write_bed(&config.after_bed, &after_site_stats).unwrap();

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

pub fn write_bed(output_filename: &String, sites: &HashSet<bamreadfilt::VCFRecord>) -> Result<(), std::io::Error>{
    use std::io::Write;
    use std::fs::File;

    let mut f = File::create(output_filename)?;

    for site in sites.iter(){
        writeln!(&mut f, "{}\t{}\t{}\t{}", site.chrm, site.pos, site.ref_b, site.alt_b)?;
    }
    Ok(())

}


