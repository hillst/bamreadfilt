extern crate bamreadfilt;
extern crate rust_htslib;


use std::env;
use std::process;
use bamreadfilt::Config;
use std::thread;
use std::time::Duration;


fn main() {
    let args: Vec<String> = env::args().collect();

    // we will need to configure our own setup here.
    let config = Config::new(&args).unwrap_or_else(|err| {
             println!("Problem parsing arguments: {}", err);
             process::exit(1);
    });

    use std::path::Path;
    use rust_htslib::bcf::Reader;
    let path = Path::new("test.vcf");
    let mut vcf_reader = Reader::from_path(path).ok().unwrap();


    for rec in vcf_reader.records(){
        let mut record = rec.ok().expect("Error reading record.");
        let allele_depth  : () = record.info(b"AD").string().unwrap().unwrap();
        println!("{:?}", allele_depth);

        //what can i do with this type?


    }


    //primary difference between this and unwrap_or_else is that since we return () in the Okay case 
    //      we have nothing we need to unwrap

    //if let Err(e) = minigrep::run(config){ // bc we return a Result type, we explicitly handle the Err case (we also could handle the () -- aka unit -- case )
    //    println!("Application error: {}", e);
    //    process::exit(1);

    //} 
}

