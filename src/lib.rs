use std::error::Error;
use std::fs::File;
use std::io::prelude::*;


/*
    usecase:
    You have a bam and you want it only to contain reads from VCF
    With no uses this is: "get all reads that contain an annotated variant in [VCF]"

    usecase advanced:
    Get all annotated variants that meet criteria:
        - bq < 20
        - mq < 20
        - mrbq < 40
        - READ ORIENTATION > 1
        - NOT IN SOFTMASK
        - fragment size ?
        - no insertions
        - no deletions ...
        - include readpairs
        etc

    programmable is better

    Functions and structs:


    Run! --
        iterate through bam
            apply filtering rules to relevant variants 


    bcf::mod::Reader::from_path

        test this on vcf and see if it works okay. if not we can stick to bcfs i guess.
 */


pub struct Config {
    pub query: String,
    pub filename: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str>{
        if args.len() < 3 {
            return Err("not enough arguments"); // why did I need return + statement here instead of just expr
        }
        let query = args[1].clone();
        let filename = args[2].clone();
        Ok(Config { query, filename })
        
    }
}

pub fn run() {
    println!("{}", "hello world");
}


mod test {
    use super::*;

    #[test]
    fn one_result() {

    }
}
