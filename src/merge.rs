use crate::*;

use sbwt::SbwtIndex;
use sbwt::SubsetSeq;
use sbwt::MergeInterleaving;

use bitvec::prelude::*;

pub fn construct_merge_plan<SS1, SS2>(sbwt_1: &SbwtIndex<SS1>, sbwt_2: &SbwtIndex<SS2>) -> MergeInterleaving
where
    SS1: SubsetSeq,
    SS2: SubsetSeq,
{
    let k = sbwt_1.k();
    let lce = LCE::new(sbwt_1, sbwt_2);

    log::info!("merging:");
    let sbwt_1_size = sbwt_1.n_sets();
    let sbwt_2_size = sbwt_2.n_sets();
    let capacity = sbwt_1_size + sbwt_2_size;

    let progress_step = capacity / 10;
    let mut checkpoint = progress_step;

    let mut s1        = BitVec::<u64, Lsb0>::with_capacity(capacity);
    let mut s2        = BitVec::<u64, Lsb0>::with_capacity(capacity);
    let mut is_dummy  = BitVec::<u64, Lsb0>::with_capacity(capacity);
    let mut is_leader = BitVec::<u64, Lsb0>::with_capacity(capacity);

    let mut previous_rank   = 0;
    let mut previous_origin = 0;
    let mut rank_1 = 1;
    let mut rank_2 = 1;

    s1.push(true);
    s2.push(true);
    is_dummy.push(true);
    is_leader.push(true);

    while rank_1 < sbwt_1_size && rank_2 < sbwt_2_size {
        let (extension, ordering) = lce.query(rank_1, 0, rank_2, 1);
        match ordering {
            std::cmp::Ordering::Less    => {
                s1.push(true);
                s2.push(false);

                is_dummy.push(lce.is_dummy(rank_1, 0));
                let (shared_suffix, _) = lce.query(rank_1, 0, previous_rank, previous_origin);
                is_leader.push(shared_suffix < k - 1);
                previous_rank   = rank_1;
                previous_origin = 0;

                rank_1 += 1;
            },
            std::cmp::Ordering::Equal   => {
                s1.push(true);
                s2.push(true);

                is_dummy.push(lce.is_dummy(rank_1, 0));
                let (shared_suffix, _) = lce.query(rank_1, 0, previous_rank, previous_origin);
                is_leader.push(shared_suffix < k - 1);
                previous_rank   = rank_1;
                previous_origin = 0;

                rank_1 += 1;
                rank_2 += 1;
            },
            std::cmp::Ordering::Greater => {
                s1.push(false);
                s2.push(true);

                is_dummy.push(lce.is_dummy(rank_2, 1));
                let (shared_suffix, _) = lce.query(rank_2, 1, previous_rank, previous_origin);
                is_leader.push(shared_suffix < k - 1);
                previous_rank   = rank_2;
                previous_origin = 1;

                rank_2 += 1;
            },
        }
        if rank_1 + rank_2 >= checkpoint {
            log::info!("done with {}/{}/{}", rank_1, rank_2, capacity);
            checkpoint += progress_step;
        }
    }
    while rank_1 < sbwt_1_size {
        s1.push(true);
        s2.push(false);
        is_dummy.push(lce.is_dummy(rank_1, 0));
        let (shared_suffix, _) = lce.query(rank_1, 0, previous_rank, previous_origin);
        is_leader.push(shared_suffix < k - 1);
        previous_rank   = rank_1;
        previous_origin = 0;
        rank_1 += 1;

        if rank_1 + rank_2 >= checkpoint {
            log::info!("done with {}/{}/{}", rank_1, rank_2, capacity);
            checkpoint += progress_step;
        }
    }

    while rank_2 < sbwt_2_size {
        s1.push(false);
        s2.push(true);
        is_dummy.push(lce.is_dummy(rank_2, 1));
        let (shared_suffix, _) = lce.query(rank_2, 1, previous_rank, previous_origin);
        is_leader.push(shared_suffix < k - 1);
        previous_rank   = rank_2;
        previous_origin = 1;
        rank_2 += 1;

        if rank_1 + rank_2 >= checkpoint {
            log::info!("done with {}/{}/{}", rank_1, rank_2, capacity);
            checkpoint += progress_step;
        }
    }

    MergeInterleaving { s1, s2, is_dummy, is_leader }
}

