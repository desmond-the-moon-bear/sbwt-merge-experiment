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
    log::info!("calculated merge:");
    let (lce, calculated_merge) = merge::construct_merge_plan(&sbwt_1, &sbwt_2);
    println!("{} {}", canonic_merge.s1.len(), calculated_merge.s1.len());
    println!("{} {}", canonic_merge.s2.len(), calculated_merge.s2.len());
    println!("{} {}", canonic_merge.is_dummy.len(), calculated_merge.is_dummy.len());
    println!("{} {}", canonic_merge.is_leader.len(), calculated_merge.is_leader.len());
    println!("{}", canonic_merge.s1 == calculated_merge.s1);
    println!("{}", canonic_merge.s2 == calculated_merge.s2);
    println!("{}", canonic_merge.is_leader == calculated_merge.is_leader);
}

fn mini_test() {
    let k = 7;
    let seq_1 = seqs![
        b"AACCCCC",
        b"AAGGGGG",
        b"TTCCCCG",
        b"GCGCGCGCG",
        b"GCGCGCA",
        b"GCGCGCT",
    ];
    let (mut sbwt_1, _) = BitPackedKmerSortingMem::new_from_vecs(&seq_1, k)
        .build_select_support(true)
        .run();

    let seq_2 = seqs![
        b"AACCCCC",
        b"CCGGGGG",
        b"CCTTTTT",
    ];
    let (mut sbwt_2, _) = BitPackedKmerSortingMem::new_from_vecs(&seq_2, k)
        .build_select_support(true)
        .run();

    let sbwts = &[&sbwt_1, &sbwt_2];

    let lce = LCE::new(&sbwt_1, &sbwt_2);

    let offset = sbwt_1.n_sets();
    let count = offset + sbwt_2.n_sets();
    for i in 1..count {
        let (index_1, origin_1) = if i < offset {
            (i, 0)
        } else {
            (i - offset, 1)
        };

        let kmer_1 = sbwts[origin_1].access_kmer(index_1);
        let kmer_1_str = str::from_utf8(&kmer_1).unwrap();
        for j in i + 1..count {
            let (index_2, origin_2) = if j < offset {
                (j, 0)
            } else {
                (j - offset, 1)
            };

            let kmer_2 = sbwts[origin_2].access_kmer(index_2);
            let kmer_2_str = str::from_utf8(&kmer_2).unwrap();

            let (extension, ordering) = lce.query(index_1, origin_1, index_2, origin_2);
            println!("{} {} {}; lce: {}", kmer_1_str, match ordering {
                std::cmp::Ordering::Less =>    "is before  ",
                std::cmp::Ordering::Equal =>   "is equal to",
                std::cmp::Ordering::Greater => "is after   ",
            }, kmer_2_str, extension);
        }
    }
}

