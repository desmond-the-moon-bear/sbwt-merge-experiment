use crate::dctree;
use crate::dctree::*;

use crate::LAQ;
use sbwt::SubsetSeq;
use sbwt::BitPackedKmerSortingMem;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

// branch_identifier:
// (type:1) [0..=k]      - branch originally connected to the root
//                         or of a loop with length greater than k
// (type:2) [k+1..=2k-1] - node in an unrolled loop
// (type:3) [2k..]       - in the subtree of a node in an unrolled loop;
//                         points to the root of the subtree

// static mut LAST_TYPE: &str = "";
// static mut LAST_RESULT: (usize, usize, usize, usize, usize) = (usize::MAX, usize::MAX, usize::MAX, usize::MAX, usize::MAX);

pub struct LongestCommonExtension {
    pub k: usize,
    pub decomp: TreeDecompositionRMQ,
    pub laq: LAQ,
}

impl LongestCommonExtension {
    pub fn new(sbwts: &[&SbwtIndexVariant]) -> LongestCommonExtension
    {
        log::info!("tree decomposition:");
        let decomp = dctree::tree_decomposition_rmq(sbwts);
        log::info!("level ancestor queries:");
        let laq = LAQ::new(&decomp.tree);

        LongestCommonExtension {
            k: sbwts[0].k(),
            laq,
            decomp,
        }
    }

    pub fn query(
        &self,
        colex_rank_1: usize,
        origin_1: usize,
        colex_rank_2: usize,
        origin_2: usize
    ) -> (usize, std::cmp::Ordering) {
        if colex_rank_1 == 0 && colex_rank_2 == 0 {
            return (self.k, std::cmp::Ordering::Equal);
        }
        let node_1 = self.decomp.mapping.node(colex_rank_1, origin_1);
        let node_2 = self.decomp.mapping.node(colex_rank_2, origin_2);
        let branch_identifier_1 = self.decomp.mapping.branch_id(colex_rank_1, origin_1);
        let branch_identifier_2 = self.decomp.mapping.branch_id(colex_rank_2, origin_2);
        self.query_nodes(node_1, branch_identifier_1, node_2, branch_identifier_2)
    }

    fn query_nodes(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        let k = self.k;

        if branch_identifier_1 <= k && branch_identifier_2 <= k {
            // Both nodes are type:1.
            return self.query_type_11(node_1, branch_identifier_1, node_2, branch_identifier_2);
        }

        // query_type_13 + query_type_12
        if branch_identifier_1 <= k || branch_identifier_2 <= k {
            // At least one node is type:1.

            let mut reverse = false;
            if branch_identifier_2 <= k {
                // Swap the nodes such that node_1 is type:1.
                reverse = true;
                std::mem::swap(&mut node_1, &mut node_2);
                std::mem::swap(&mut branch_identifier_1, &mut branch_identifier_2);
            }

            // node_1 is type:1, node_2 is type:2 or type:3.

            let mut subtree_extension = 0;
            if branch_identifier_2 >= 2 * k {
                // node_2 is type:3.
                subtree_extension += self.rm_query(node_1, node_2);
                if subtree_extension >= k {
                    return (k, std::cmp::Ordering::Equal);
                }
                if subtree_extension < self.decomp.tree[node_2].depth {
                    // We have not reached the loop.
                    let mut ordering = if subtree_extension == branch_identifier_1 {
                        // node_1 corresponds to a dummy, therefore, it is smaller since node_2 is
                        // guaranteed to be a non-dummy.
                        std::cmp::Ordering::Less
                    } else {
                        // Neither node is a dummy and the extension is shorter k, therefore, we
                        // compare the ancestors.
                        self.compare_ancestors(node_1, node_2, subtree_extension)
                    };
                    if reverse {
                        ordering = ordering.reverse();
                    }
                    return (subtree_extension, ordering);
                }
                // Move node_1 to the corresponding ancestor.
                node_1 = self.la_query(node_1, subtree_extension);
                branch_identifier_1 = self.find_branch_identifier(node_1);
                // Move node_2 to the corresponding type:2 node.
                node_2 = branch_identifier_2 - 2 * k;
                branch_identifier_2 = self.find_branch_identifier(node_2);
            }

            // node_2 is type:2.
            debug_assert!(branch_identifier_1 <= k);
            debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);
            let first_loop_extension = self.rm_query(node_1, node_2);
            if first_loop_extension + subtree_extension >= k {
                return (k, std::cmp::Ordering::Equal);
            }
            if first_loop_extension < self.decomp.tree[node_2].depth {
                let mut ordering = if first_loop_extension == branch_identifier_1 {
                    std::cmp::Ordering::Less
                } else {
                    self.compare_ancestors(node_1, node_2, first_loop_extension)
                };
                if reverse {
                    ordering = ordering.reverse();
                }
                return (first_loop_extension + subtree_extension, ordering);
            }

            let loop_2_length = self.decomp.tree[node_2].depth + k - branch_identifier_2;
            let depth_in_loop_2 = self.decomp.tree[node_2].depth - loop_2_length;
            debug_assert!(loop_2_length + depth_in_loop_2 == self.decomp.tree[node_2].depth);

            // <---------> loop length
            //             <------> depth in loop
            // [] .. .. [] [] .. [] ..
            //                   ^^ - node_2
            //
            // --------------------------------------------------------------------
            //
            //                <---------> another possible instance of the loop
            //                            <---------> loop length
            //                                        <------> depth in loop
            // .. .. .. .. {} [] .. .. [] [] .. .. [] [] .. []
            //             ^^          ^^          ^^       ^^ - node_1
            //             ^^          ^^          ^^ - sample ancestor
            //             ^^          ^^ - shifted ancestor
            //             ^^ - different character
            // 
            // shift ancestor < sample ancestor <=> node_1 < node_2
            // shift ancestor > sample ancestor <=> node_1 > node_2
            //

            let sample_ancestor  = self.la_query(         node_1, depth_in_loop_2);
            let shifted_ancestor = self.la_query(sample_ancestor, loop_2_length);
            let sample_branch_identifier  = self.find_branch_identifier(sample_ancestor);
            let shifted_branch_identifier = self.find_branch_identifier(shifted_ancestor);
            let (extension, mut ordering) = self.query_type_11(
                shifted_ancestor,
                shifted_branch_identifier,
                sample_ancestor,
                sample_branch_identifier,
            );

            let total_extension = extension + self.decomp.tree[node_2].depth + subtree_extension;
            if total_extension >= k {
                return (k, std::cmp::Ordering::Equal);
            }
            if reverse {
                ordering = ordering.reverse();
            }
            return (total_extension, ordering);
        }

        if 2 * k <= branch_identifier_1 && 2 * k <= branch_identifier_2 {
            // Both nodes are type:3.
            return self.query_type_33(node_1, branch_identifier_1, node_2, branch_identifier_2);
        }

        self.dispatch_query_23_or_22(node_1, branch_identifier_1, node_2, branch_identifier_2)
    }

    fn query_type_11(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        // unsafe {
        //     LAST_TYPE = "11";
        // }
        let k = self.k;
        let extension = self.rm_query(node_1, node_2);
        if extension >= k {
            return (k, std::cmp::Ordering::Equal);
        }
        if extension == branch_identifier_1 || extension == branch_identifier_2 {
            // At least one is a dummy k-mer since [extension < k].
            if branch_identifier_1 == branch_identifier_2 {
                // Both are dummy k-mers.
                return (extension, std::cmp::Ordering::Equal);
            }
            // Exactly one is a dummy k-mer. The dummy k-mer comes before the non-dummy k-mer.
            return (extension, branch_identifier_1.cmp(&branch_identifier_2));
        }
        let ordering = self.compare_ancestors(node_1, node_2, extension);
        (extension, ordering)
    }

    fn query_type_33(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        // unsafe {
        //     LAST_TYPE = "33";
        // }
        let k = self.k;

        debug_assert!(2 * k <= branch_identifier_1);
        debug_assert!(2 * k <= branch_identifier_2);

        let subtree_extension = self.rm_query(node_1, node_2);

        if subtree_extension < self.decomp.tree[node_1].depth && subtree_extension < self.decomp.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, subtree_extension);
            return (subtree_extension, ordering);
        }

        if subtree_extension >= self.decomp.tree[node_1].depth {
            node_1 = branch_identifier_1 - 2 * k;
            branch_identifier_1 = self.find_branch_identifier(node_1);
            debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        } else {
            node_1 = self.la_query(node_1, subtree_extension);
            // branch_identifier_1 = self.find_branch_identifier(node_1);
            debug_assert!(2 * k <= branch_identifier_1);
        }

        if subtree_extension >= self.decomp.tree[node_2].depth {
            node_2 = branch_identifier_2 - 2 * k;
            branch_identifier_2 = self.find_branch_identifier(node_2);
            debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);
        } else {
            node_2 = self.la_query(node_2, subtree_extension);
            // branch_identifier_2 = self.find_branch_identifier(node_2);
            debug_assert!(2 * k <= branch_identifier_2);
        }

        let (extension, ordering) = self.dispatch_query_23_or_22(
            node_1,
            branch_identifier_1,
            node_2,
            branch_identifier_2
        );

        if extension + subtree_extension >= k {
            return (k, std::cmp::Ordering::Equal);
        }

        (extension + subtree_extension, ordering)
    }

    fn dispatch_query_23_or_22(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        let k = self.k;
        if branch_identifier_1 >= 2 * k || branch_identifier_2 >= 2 * k {
            // One node is type:3 and one node is type:2.
            let mut reverse = false;
            if branch_identifier_1 >= 2 * k {
                // node_1 is type:3. We want node_2 to be type 3.
                reverse = true;
                std::mem::swap(&mut node_1, &mut node_2);
                std::mem::swap(&mut branch_identifier_1, &mut branch_identifier_2);
            }
            let (extension, mut ordering) = self.query_type_23(
                node_1,
                branch_identifier_1,
                node_2,
                branch_identifier_2
            );
            if reverse {
                ordering = ordering.reverse();
            }
            return (extension, ordering);
        }

        // Both nodes are type:2.
        self.query_type_22(
            node_1,
            branch_identifier_1,
            node_2,
            branch_identifier_2
        )
    }

    fn query_type_23(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        // unsafe {
        //     LAST_TYPE = "23"
        // }
        let k = self.k;

        debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        debug_assert!(2 * k <= branch_identifier_2);

        let mut subtree_extension = self.rm_query(node_1, node_2);

        if subtree_extension < self.decomp.tree[node_1].depth && subtree_extension < self.decomp.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, subtree_extension);
            return (subtree_extension, ordering);
        }

        if subtree_extension < self.decomp.tree[node_1].depth || subtree_extension >= self.decomp.tree[node_2].depth {
            // (subtree_extension < self.rmq.tree[node_1].depth and subtree_extension >= self.rmq.tree[node_2].depth)
            // or (subtree_extension >= self.rmq.tree[node_1].depth and subtree_extension >= self.rmq.tree[node_2].depth)
            let jump_distance = if subtree_extension < self.decomp.tree[node_1].depth {
                subtree_extension
            } else {
                let loop_1_length = self.decomp.tree[node_1].depth + k - branch_identifier_1;
                self.decomp.tree[node_1].depth - loop_1_length // depth_in_loop_1
            };
            let node_1_jump = self.la_query(node_1, jump_distance);
            let (node_1, branch_identifier_1) = self.find_canonic_node_and_branch_identifier(node_1_jump);
            debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);

            node_2 = branch_identifier_2 - 2 * k;
            branch_identifier_2 = self.find_branch_identifier(node_2);
            debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);

            let (extension, ordering) = self.query_type_22(node_1, branch_identifier_1, node_2, branch_identifier_2);
            if subtree_extension + extension >= k {
                return (k, std::cmp::Ordering::Equal);
            }
            return (subtree_extension + extension, ordering);
        }

        // subtree_extension >= self.rmq.tree[node_1].depth and subtree_extension < self.rmq.tree[node_2].depth

        let loop_1_length = self.decomp.tree[node_1].depth + k - branch_identifier_1;
        let depth_in_loop_1 = self.decomp.tree[node_1].depth - loop_1_length;

        debug_assert!(loop_1_length + depth_in_loop_1 == self.decomp.tree[node_1].depth);
        debug_assert!(self.decomp.tree[node_1].depth < self.decomp.tree[node_2].depth);

        let sample_ancestor  = self.la_query(         node_2, depth_in_loop_1);
        let shifted_ancestor = self.la_query(sample_ancestor, loop_1_length);
        let shift_extension = self.rm_query(sample_ancestor, shifted_ancestor);

        let total_subtree_extension = shift_extension + self.decomp.tree[node_1].depth;
        if total_subtree_extension < self.decomp.tree[node_2].depth {
            // non-loop node === shifted ancestor
            let ordering = self.compare_ancestors(sample_ancestor, shifted_ancestor, shift_extension);
            return (total_subtree_extension, ordering);
        }

        let node_1_jump = self.la_query(node_1, depth_in_loop_1);
        let (node_1, branch_identifier_1) = self.find_canonic_node_and_branch_identifier(node_1_jump);
        debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);

        node_2 = branch_identifier_2 - 2 * k;
        branch_identifier_2 = self.find_branch_identifier(node_2);
        debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);

        let (loop_extension, ordering) = self.query_type_22(node_1, branch_identifier_1, node_2, branch_identifier_2);
        debug_assert!(loop_extension + total_subtree_extension < k);

        (loop_extension + total_subtree_extension, ordering)
    }

    fn query_type_22(
        &self,
        mut node_1: usize,
        mut branch_identifier_1: usize,
        mut node_2: usize,
        mut branch_identifier_2: usize
    ) -> (usize, std::cmp::Ordering) {
        // unsafe {
        //     LAST_TYPE = "23"
        // }
        let k = self.k;

        debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);

        let mut loop_1_length   = self.decomp.tree[node_1].depth + k - branch_identifier_1;
        let mut depth_in_loop_1 = self.decomp.tree[node_1].depth - loop_1_length;
        let mut loop_2_length   = self.decomp.tree[node_2].depth + k - branch_identifier_2;
        let mut depth_in_loop_2 = self.decomp.tree[node_2].depth - loop_2_length;

        let initial_extension = self.rm_query(node_1, node_2);

        if initial_extension < self.decomp.tree[node_1].depth && initial_extension < self.decomp.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, initial_extension);
            return (initial_extension, ordering);
        }

        if initial_extension >= self.decomp.tree[node_1].depth && initial_extension >= self.decomp.tree[node_2].depth {
            return (k, std::cmp::Ordering::Equal);
        }

        let mut reverse = false;
        if initial_extension >= self.decomp.tree[node_2].depth {
            reverse = true;
            std::mem::swap(&mut node_1,              &mut node_2);
            std::mem::swap(&mut branch_identifier_1, &mut branch_identifier_2);
            std::mem::swap(&mut loop_1_length,       &mut loop_2_length);
            std::mem::swap(&mut depth_in_loop_1,     &mut depth_in_loop_2);
        }

        // (initial_extension >= self.rmq.tree[node_1].depth and initial_extension <  self.rmq.tree[node_2].depth)

        let sample_ancestor  = self.find_canonic_node(self.la_query(         node_2, depth_in_loop_1));
        let shifted_ancestor = self.find_canonic_node(self.la_query(sample_ancestor, loop_1_length));
        let shift_extension = self.rm_query(sample_ancestor, shifted_ancestor);
        
        // longer loop === shifted_ancestor
        let mut ordering = self.compare_ancestors(sample_ancestor, shifted_ancestor, shift_extension);
        if reverse {
            ordering = ordering.reverse();
        }

        (shift_extension + self.decomp.tree[node_1].depth, ordering)
    }

    #[inline]
    fn la_query(&self, node: usize, climb: usize) -> usize {
        self.laq.query(&self.decomp.tree, node, climb)
    }

    #[inline]
    fn rm_query(&self, node_1: usize, node_2: usize) -> usize {
        let name_1 = self.decomp.tree[node_1].name;
        let name_2 = self.decomp.tree[node_2].name;
        if name_1 == name_2 {
            return self.decomp.rmq.name_lengths()[name_1]
        }
        let smaller_name = name_1.min(name_2);
        let bigger_name = name_1 + name_2 - smaller_name;
        self.decomp.rmq.query(smaller_name + 1, bigger_name)
    }

    #[inline]
    fn compare_ancestors(
        &self,
        node_1: usize,
        node_2: usize,
        extension: usize,
    ) -> std::cmp::Ordering {
        let parent_1 = self.la_query(node_1, extension);
        let parent_2 = self.la_query(node_2, extension);
        let parent_1_name = self.decomp.tree[parent_1].name;
        let parent_2_name = self.decomp.tree[parent_2].name;
        // unsafe {
        //     LAST_RESULT = (parent_1, parent_2, parent_1_name, parent_2_name, extension)
        // }
        parent_1_name.cmp(&parent_2_name)
    }

    #[inline]
    fn find_branch_identifier(&self, node: usize) -> usize {
        let (colex_rank, origin) = self.decomp.mapping.id_to_colex_rank(self.decomp.tree[node].id);
        self.decomp.mapping.branch_id(colex_rank, origin)
    }

    #[inline]
    fn find_canonic_node(&self, node: usize) -> usize {
        self.decomp.mapping.id_to_canonic_node(self.decomp.tree[node].id)
    }

    #[inline]
    fn find_canonic_node_and_branch_identifier(&self, node: usize) -> (usize, usize) {
        let id = self.decomp.tree[node].id;
        let (colex_rank, origin) = self.decomp.mapping.id_to_colex_rank(id);
        self.decomp.mapping.node_and_branch_id(colex_rank, origin)
    }

    pub fn is_dummy(
        &self,
        colex_rank: usize,
        origin: usize,
    ) -> bool {
        self.decomp.mapping.branch_id(colex_rank, origin) < self.k
    }
}

pub fn test_lce_from_sbwts(sbwts: &[&SbwtIndexVariant]) {
    let k = sbwts[0].k();
    let lce = LongestCommonExtension::new(sbwts);

    log::info!("testing laq:");
    fn level_ancestor_query(tree: &[Node], mut node: usize, mut climb: usize) -> usize {
        while climb > 0 {
            node = tree[node].parent;
            climb -= 1;
        }
        node
    }

    for (node_index, node) in lce.decomp.tree.iter().enumerate() {
        // let max_depth = node.depth;
        let max_depth = (2 * k).min(lce.decomp.tree[node_index].depth);
        for climb in 0..max_depth {
            let invariant = level_ancestor_query(&lce.decomp.tree, node_index, climb);
            let found     = lce.laq.query       (&lce.decomp.tree, node_index, climb);
            assert_eq!(invariant, found, "node_index: {}; climb: {}", node_index, climb);
        }
    }

    log::info!("testing lce:");
    let count = lce.decomp.mapping.sbwt_to_tree.len();

    let mut kmer_i_buf = Vec::<u8>::new();
    let mut kmer_j_buf = Vec::<u8>::new();
    let mut kmer_i_parent_buf = Vec::<u8>::new();
    let mut kmer_j_parent_buf = Vec::<u8>::new();

    let mut origin_i = 0;
    for i in 0..count {
        if i >= lce.decomp.mapping.sbwt_offsets[origin_i + 1] {
            origin_i += 1;
        }
        let index_i = i - lce.decomp.mapping.sbwt_offsets[origin_i];

        kmer_i_buf.clear();
        sbwts[origin_i].push_kmer_to_vec(index_i, &mut kmer_i_buf);
        let kmer_i_str = str::from_utf8(&kmer_i_buf).unwrap();

        let mut origin_j = 0;
        for j in 0..count {
            if j >= lce.decomp.mapping.sbwt_offsets[origin_j + 1] {
                origin_j += 1;
            }
            let index_j = j - lce.decomp.mapping.sbwt_offsets[origin_j];

            if j == i {
                continue;
            }

            kmer_j_buf.clear();
            sbwts[origin_j].push_kmer_to_vec(index_j, &mut kmer_j_buf);
            let kmer_j_str = str::from_utf8(&kmer_j_buf).unwrap();

            let (extension, ordering) = lce.query(index_i, origin_i, index_j, origin_j);
            let mut true_extension = 0;
            for j in (0..k).rev() {
                if kmer_i_buf[j] != kmer_j_buf[j] {
                    break;
                }
                true_extension += 1;
            }
            let correct = if ordering == std::cmp::Ordering::Equal {
                kmer_i_buf == kmer_j_buf
            } else {
                true_extension == extension &&
                kmer_i_buf[k-extension..k] ==  kmer_j_buf[k-extension..k] &&
                kmer_i_buf[k-extension-1].cmp(&kmer_j_buf[k-extension-1]) == ordering
            };
            let (node_i, branch_identifier_i) = lce.decomp.mapping.node_and_branch_id(index_i, origin_i);
            let (node_j, branch_identifier_j) = lce.decomp.mapping.node_and_branch_id(index_j, origin_j);
            if !correct {
                // unsafe {
                //     let thing = LAST_TYPE;
                //     println!("{}", thing);
                //     let thing = LAST_RESULT;
                //     println!("{:?}", thing);
                // }
                println!("index_i: {}, origin_i: {}, node_i: {}", index_i, origin_i, node_i);
                println!("index_j: {}, origin_j: {}, node_j: {}", index_j, origin_j, node_j);
                for climb in 1..=extension {
                    kmer_i_parent_buf.clear();
                    kmer_j_parent_buf.clear();
                    let node_i_parent = lce.decomp.tree[lce.la_query(node_i, climb)].id;
                    let node_j_parent = lce.decomp.tree[lce.la_query(node_j, climb)].id;
                    let (node_i_rank, parent_origin_i) = lce.decomp.mapping.id_to_colex_rank(node_i_parent);
                    let (node_j_rank, parent_origin_j) = lce.decomp.mapping.id_to_colex_rank(node_j_parent);
                    sbwts[parent_origin_i].push_kmer_to_vec(node_i_rank, &mut kmer_i_parent_buf);
                    sbwts[parent_origin_j].push_kmer_to_vec(node_j_rank, &mut kmer_j_parent_buf);
                    println!("{}: {} {}",
                        climb,
                        str::from_utf8(&kmer_i_parent_buf).unwrap(),
                        str::from_utf8(&kmer_j_parent_buf).unwrap(),
                    );
                }
            }
            assert!(
                correct,
                "{} {} {} {:?} {} {}", kmer_i_str, kmer_j_str, extension, ordering,
                branch_identifier_i, branch_identifier_j,
            );
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use super::LongestCommonExtension as LCE;

    fn test_lce_from_seqs(k: usize, seq_1: Vec<Vec<u8>>, seq_2: Vec<Vec<u8>>) 
    {
        let (mut sbwt_1, _) = BitPackedKmerSortingMem::new_from_vecs(&seq_1, k)
            .build_select_support(true)
            .run();
        let (mut sbwt_2, _) = BitPackedKmerSortingMem::new_from_vecs(&seq_2, k)
            .build_select_support(true)
            .run();

        let sbwts = &[
            &SbwtIndexVariant::SubsetMatrix(sbwt_1),
            &SbwtIndexVariant::SubsetMatrix(sbwt_2),
        ];

        test_lce_from_sbwts(sbwts);
    }

    #[test]
    fn lce_1() {
        let k = 7;
        let seq_1 = seqs![
            b"AACCCCC",
            b"AAGGGGG",
            b"TTCCCCG",
            b"GCGCGCGCG",
            b"GCGCGCA",
            b"GCGCGCT",
        ];
        let seq_2 = seqs![
            b"AACCCCC",
            b"CCGGGGG",
            b"CCTTTTT",
            b"ACACCACCCACCCC",
            b"ACACCAGGGAGGGG",
            b"ACAGGAGGGAGGGG",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_2() {
        let k = 7;
        let seq_1 = seqs![
            b"ATATATATA",
            b"ATATACGCG",
        ];
        let seq_2 = seqs![
            b"CGCGCGCGCG",
            b"CGCGCATAT",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_3() {
        let k = 7;
        let seq_1 = seqs![
            b"ATATATATATCGCGCG",
        ];
        let seq_2 = seqs![
            b"CGCGCGCGCGATATAT",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_4() {
        let k = 7;
        let seq_1 = seqs![
            b"GATATATGATATAT",
        ];
        let seq_2 = seqs![
            b"ATATATAT",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_5() {
        let k = 7;
        let seq_1 = seqs![
            b"CATCATCATCATCAT",
        ];
        let seq_2 = seqs![
            b"CTACTACTACTACTA",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_6() {
        let k = 7;
        let seq_1 = seqs![
            b"CATCATCATCATCAT",
        ];
        let seq_2 = seqs![
            b"CATCATCATCATCAT",
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }

    #[test]
    fn lce_7() {
        let k = 7;
        let seq_1 = seqs![
            b"GACACGACACGACACGACAC",
            b"GTGTGTGTGTGTGT"
        ];
        let seq_2 = seqs![
            b"ACACACACACACAC",
            b"CGTGTCGTGTCGTGTCGTGT"
        ];
        test_lce_from_seqs(k, seq_1, seq_2);
    }
}
