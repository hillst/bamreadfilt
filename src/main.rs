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
        mscp_run(config, 2 as usize);
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
    use std::io::Write;
    use std::fs::File;
    
    let bam_filename = config.bam_filename.clone();
    let vcf_filename= config.vcf_filename.clone();
    
    //TODO duplicated code
    let vcf_path = Path::new(&config.vcf_filename);
    let mut vcf_itr = Reader::from_path(vcf_path).expect("cannot open vcf"); 

    let path = Path::new(&bam_filename);//why is it only using move semantics lol
    let bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");
    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;
   
    let mut out = match config.use_stdout {
        false => bam::Writer::from_path("/dev/null", &head).expect("Error opening bam for writing."), 
        true => bam::Writer::from_stdout(&head).expect("Error opening stdout for writing."),
    };

    let mut raw_stats_fd = File::create(config.raw_stats_fn.clone()).unwrap();

    
    let mut before_stats = bamreadfilt::SiteStats::new();
    let mut after_stats = bamreadfilt::SiteStats::new();

    let mut before_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();
    let mut after_site_stats: HashSet<bamreadfilt::VCFRecord> = HashSet::new();

    let mut vec_list: Vec<Vec<rust_htslib::bcf::Record>> = vec![];


    let mut full_list: Vec<rust_htslib::bcf::Record> = vec![];
    let mut master_list: Vec<rust_htslib::bcf::Record> = vec![]; // includes downsampling

    if config.max_sites > 0 && config.max_sites < n_threads as i64 {
        eprintln!("[warn]: Downsampled sites must be larger than the number of threads.");
        return;
    }
    
    for (i, entry) in vcf_itr.records().enumerate() {
        // I'm almost certain we will want these to be continuous in memory or atleast adjacent
        master_list.push(entry.ok().unwrap());
    }
    
    // if master list is longer than n_threads we ok
    //   otherwise decrease the number of threads? Or???? panic? warn the user?
    if master_list.len() < n_threads {
        eprintln!("[warn]: Empty VCF or number of mutations < number of threads. Are you using downsampling? You MUST have more sites than threads.");
        return;
    }

    if config.max_sites > 0 && config.max_sites < master_list.len() as i64 {
        eprintln!("Max Sites: {}  Total Sites: {}", config.max_sites, master_list.len() as i64);
        //master_list = master_list.head(config.max_sites as usize);
        master_list.truncate(config.max_sites as usize);
    }
        

    
    for _ in 0..n_threads{
        vec_list.push(vec![]);
    }

    for (i, entry) in master_list.into_iter().enumerate() {
        // I'm almost certain we will want these to be continuous in memory or atleast adjacent
        vec_list[i % n_threads].push(entry); //wait i wanna move this
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
        let _vp = Path::new(&vcf_filename); 
        let _vr = rust_htslib::bcf::Reader::from_path(_vp).ok().expect("Error opening vcf.");
        let vcf_head: rust_htslib::bcf::header::HeaderView = _vr.header().clone(); //jfc
        for (i, (mut item, _vcf, _pir)) in rx.iter()
            .filter(| &(ref i, v, p)| {
               i.seq()[p] == v.alt_b 
            })
            .map(   | (_i, _v, p) | { 
                before_stats.collect_stats_stream(&_i, &_v, p); //NOTE this works because it is the consumer, so we have no race conditions at this point
                before_site_stats.insert(_v.clone());
                (_i, _v, p)
            })
            .filter(| &(ref _i, _v, p)| { ! config.remove_duplicates || ! _i.is_duplicate() }) //rmdup
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
                /* Prototype:
                 *    This approach woudl let us offload the labelling (tumor-normal) to the end user for our usecase
                        it also would allow more general usage outside of our own.
                 * 
                *  if config.label_tag OR for tag in config.tags (????)

                    writeline(...., tag)
                */

                item.push_aux( "B0".as_bytes(), &Aux::Integer(_vcf.pos as i64) );
                item.push_aux( "B1".as_bytes(), &Aux::Char(_vcf.ref_b) );
                item.push_aux( "B2".as_bytes(), &Aux::Char(_vcf.alt_b) );
                item.push_aux( "B3".as_bytes(), &Aux::Integer(_pir as i64) );

            
                if config.raw_stats_fn != "/dev/null".to_string() { 
                    //TODO write a test here, _vcf.pos alone was giving us incorrect results :)
                    if config.tag != "".to_string(){
                        
                        let test: &[u8] = config.tag.as_bytes();
                        let aux_val = match item.aux(test){
                            Some(v) => v,
                            None => {continue;}
                        };
                        let another_test = aux_val.string();
                    //exact same writeline as before but with the specified tag included as well
                        writeln!(raw_stats_fd,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",  
                                     std::str::from_utf8(vcf_head.rid2name(_vcf.chrm)).unwrap(), _vcf.pos , _vcf.ref_b as char, _vcf.alt_b as char,  // header info
                                     std::str::from_utf8( &item.qname()).unwrap(),                 // 
                                     bamreadfilt::mrbq(&item), &item.qual()[_pir] , item.mapq(), _pir, item.insert_size() , item.is_reverse() as i32,
                                     item.is_first_in_template() as i32,
                                     std::str::from_utf8(another_test).unwrap(), 
                        );
                    }
                    else {
                        writeln!(raw_stats_fd,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",  
                                     std::str::from_utf8(vcf_head.rid2name(_vcf.chrm)).unwrap(), _vcf.pos , _vcf.ref_b as char, _vcf.alt_b as char,  // header info
                                     std::str::from_utf8( &item.qname()).unwrap(),                 // 
                                     bamreadfilt::mrbq(&item), &item.qual()[_pir] , item.mapq(), _pir, item.insert_size(), item.is_reverse() as i32,
                                     item.is_first_in_template() as i32,
                        );
                    }
                }
                else {
                    //if we arent instats mode, we write a bam!
                    out.write(&item).unwrap(); // writes filtered reads.
                }
            }
            

        eprintln!("Done.");
        //write_statistics(&config.stats, before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
        //let _vp = Path::new(&vcf_filename); 
        //let _vr = rust_htslib::bcf::Reader::from_path(_vp).ok().expect("Error opening vcf.");
        //let vcf_head: rust_htslib::bcf::header::HeaderView = _vr.header().clone(); //jfc
        if config.stats_before{
            eprintln!("WARNING: using new version of brf which does not record stats by default!");
            write_bed(&config.before_bed, &before_site_stats, &vcf_head).unwrap();  
        }
        if config.stats_after{
            eprintln!("WARNING: using new version of brf which does not record stats by default!");
            write_bed(&config.after_bed, &after_site_stats, &vcf_head).unwrap();
        }

    });
    eprintln!("Cleaning up.");

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
    let vcf_filename= config.vcf_filename.clone();

    let _vp = Path::new(&vcf_filename); 
    let _vr = rust_htslib::bcf::Reader::from_path(_vp).ok().expect("Error opening vcf.");
    let vcf_head: rust_htslib::bcf::header::HeaderView = _vr.header().clone(); //jfc
    //only used for building header
    let path = Path::new(&bam_filename); //why is it only using move semantics lol
    let bam = bam::IndexedReader::from_path(path).ok().expect("Error opening bam.");
    let head = bam::Header::from_template(bam.header()); // wonder if i gotta clone dis;

    //sometimes we dont want stdout. what do?
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
                                        .filter(| &(ref _i, _v, p)| { ! config.remove_duplicates || ! _i.is_duplicate() }) //rmdup
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
        item.push_aux( "B0".as_bytes(), &Aux::Integer(_vcf.pos as i64) );
        item.push_aux( "B1".as_bytes(), &Aux::Char(_vcf.ref_b) );
        item.push_aux( "B2".as_bytes(), &Aux::Char(_vcf.alt_b) );
        item.push_aux( "B3".as_bytes(), &Aux::Integer(_pir as i64) );

        if config.use_stdout {
            out.write(&item).unwrap(); // writes filtered reads.
        }
    }

    write_statistics(&config.stats.to_string(), before_stats, after_stats, before_site_stats.len() as u64, after_site_stats.len() as u64).unwrap();
    write_bed(&config.before_bed, &before_site_stats, &vcf_head).unwrap();
    write_bed(&config.after_bed, &after_site_stats, &vcf_head).unwrap();

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

pub fn write_bed(output_filename: &String, sites: &HashSet<bamreadfilt::VCFRecord>, vcf_head: &rust_htslib::bcf::header::HeaderView) -> Result<(), std::io::Error>{
    use std::io::Write;
    use std::fs::File;

    let mut f = File::create(output_filename)?;

    for site in sites.iter(){
        writeln!(&mut f, "{}\t{}\t{}\t{}", std::str::from_utf8(vcf_head.rid2name(site.chrm)).unwrap(), site.pos, site.ref_b as char , site.alt_b as char)?;
    }
    Ok(())

}


