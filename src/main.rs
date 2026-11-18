#![allow(unused)]

macro_rules! seqs {
    ($($seq:expr),* $(,)?) => {
        vec![$($seq.to_vec(),)*]
    };
}

pub mod rmq;
pub mod level;
pub mod lce;
pub mod dctree;
pub mod merge;

use rmq::RangeMinimumQuery as RMQ;
use level::LevelAncestorSupport as LAQ;
use lce::LongestCommonExtension as LCE;

use sbwt::BitPackedKmerSortingMem;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

const PATH1: &str = "./data/v1.sbwt";
const PATH2: &str = "./data/v2.sbwt";

fn load_index(path: &str) -> SbwtIndexVariant {
    let mut reader = std::io::BufReader::new(std::fs::File::open(path).unwrap());
    SbwtIndexVariant::load(&mut reader).unwrap()
}

fn main() {
    env_logger::init();
    merge_data_test();
}

fn merge_data_test() {
    let mut sbwt_1 = load_index(PATH1);
    let mut sbwt_2 = load_index(PATH2);
    sbwt_1.build_select();
    sbwt_2.build_select();

    let sbwts = &[&sbwt_1, &sbwt_2];

    log::info!("calculated merge:");
    let (lce, calculated_merge) = merge::construct_merge_plan(sbwts);

    let sbwt_1 = match sbwt_1 {
        SbwtIndexVariant::SubsetMatrix(sbwt_index) => sbwt_index,
        SbwtIndexVariant::SubsetCorrectionSets(sbwt_index) => unreachable!(),
    };
    let sbwt_2 = match sbwt_2 {
        SbwtIndexVariant::SubsetMatrix(sbwt_index) => sbwt_index,
        SbwtIndexVariant::SubsetCorrectionSets(sbwt_index) => unreachable!(),
    };
    log::info!("canonic merge:");
    let canonic_merge = sbwt::MergeInterleaving::new(&sbwt_1, &sbwt_2, false, 7);

    println!("{}", canonic_merge.s1 == calculated_merge.s1);
    println!("{}", canonic_merge.s2 == calculated_merge.s2);
    println!("{}", canonic_merge.is_dummy == calculated_merge.is_dummy);
    println!("{}", canonic_merge.is_leader == calculated_merge.is_leader);
}

