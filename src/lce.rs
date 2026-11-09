use super::*;

// branch_identifier:
// (type:1) [0..=k]      - branch originally connected to the root
//                         or of a loop with length greater than k
// (type:2) [k+1..=2k-1] - node in an unrolled loop
// (type:3) [2k..]       - in the subtree of a node in an unrolled loop;
//                         points to the root of the subtree

pub struct LongestCommonExtension {
    k: usize,
    offset: usize,
    rmq: TreeDecompositionRMQ,
    laq: LAQ,
}

impl LongestCommonExtension {
    pub fn new<SS1, SS2>(sbwt_1: &SbwtIndex<SS1>, sbwt_2: &SbwtIndex<SS2>) -> LongestCommonExtension
    where
        SS1: SubsetSeq,
        SS2: SubsetSeq,
    {
        let rmq = range_minimum_query(sbwt_1, sbwt_2);
        let laq = LAQ::new(&rmq.tree);

        LongestCommonExtension {
            k: sbwt_1.k(),
            offset: sbwt_1.n_sets(),
            laq,
            rmq,
        }
    }

    pub fn query(
        &self,
        sbwt_index_1: usize,
        sbwt_origin_1: usize,
        sbwt_index_2: usize,
        sbwt_origin_2: usize
    ) -> (usize, std::cmp::Ordering) {
        let mapping_1 = if sbwt_origin_1 == 0 { &self.rmq.sbwt_1_mapping } else { &self.rmq.sbwt_2_mapping };
        let mapping_2 = if sbwt_origin_2 == 0 { &self.rmq.sbwt_1_mapping } else { &self.rmq.sbwt_2_mapping };
        let node_1 = mapping_1.sbwt_to_tree[sbwt_index_1];
        let node_2 = mapping_2.sbwt_to_tree[sbwt_index_2];
        let branch_identifier_1 = mapping_1.branch_identifier[sbwt_index_1];
        let branch_identifier_2 = mapping_2.branch_identifier[sbwt_index_2];
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

        let mut reverse = false;
        // query_type_13 + query_type_12
        if branch_identifier_1 <= k || branch_identifier_2 <= k {
            // At least one node is type:1.

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
                let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
                subtree_extension += self.rmq.rmq.query(rmq_start, rmq_end);
                if subtree_extension >= k {
                    return (k, std::cmp::Ordering::Equal);
                }
                if subtree_extension < self.rmq.tree[node_2].depth {
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
            debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);
            debug_assert!(branch_identifier_1 <= k);
            let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
            let first_loop_extension = self.rmq.rmq.query(rmq_start, rmq_end);
            if first_loop_extension + subtree_extension >= k {
                return (k, std::cmp::Ordering::Equal);
            }
            if first_loop_extension < self.rmq.tree[node_2].depth {
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

            let loop_2_length = self.rmq.tree[node_2].depth + k - branch_identifier_2;
            let depth_in_loop_2 = self.rmq.tree[node_2].depth - loop_2_length;
            debug_assert!(loop_2_length + depth_in_loop_2 == self.rmq.tree[node_2].depth);

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

            let total_extension = extension + depth_in_loop_2 + subtree_extension;
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
        let k = self.k;
        let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
        let extension = self.rmq.rmq.query(rmq_start, rmq_end);
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
        let k = self.k;

        debug_assert!(2 * k <= branch_identifier_1);
        debug_assert!(2 * k <= branch_identifier_2);

        let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
        let subtree_extension = self.rmq.rmq.query(rmq_start, rmq_end);

        if subtree_extension < self.rmq.tree[node_1].depth && subtree_extension < self.rmq.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, subtree_extension);
            return (subtree_extension, ordering);
        }

        if subtree_extension >= self.rmq.tree[node_1].depth {
            node_1 = branch_identifier_1 - 2 * k;
            branch_identifier_1 = self.find_branch_identifier(node_1);
            debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        } else {
            node_1 = self.la_query(node_1, subtree_extension);
            // branch_identifier_1 = self.find_branch_identifier(node_1);
            debug_assert!(2 * k <= branch_identifier_1);
        }

        if subtree_extension >= self.rmq.tree[node_2].depth {
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
        let k = self.k;

        debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        debug_assert!(2 * k <= branch_identifier_2);

        let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
        let mut subtree_extension = self.rmq.rmq.query(rmq_start, rmq_end);

        if subtree_extension < self.rmq.tree[node_1].depth && subtree_extension < self.rmq.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, subtree_extension);
            return (subtree_extension, ordering);
        }

        if subtree_extension < self.rmq.tree[node_1].depth || subtree_extension >= self.rmq.tree[node_2].depth {
            // (subtree_extension < self.rmq.tree[node_1].depth and subtree_extension >= self.rmq.tree[node_2].depth)
            // or (subtree_extension >= self.rmq.tree[node_1].depth and subtree_extension >= self.rmq.tree[node_2].depth)
            let jump_distance = if subtree_extension < self.rmq.tree[node_1].depth {
                subtree_extension
            } else {
                let loop_1_length = self.rmq.tree[node_1].depth + k - branch_identifier_1;
                self.rmq.tree[node_1].depth - loop_1_length // depth_in_loop_1
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

        let loop_1_length = self.rmq.tree[node_1].depth + k - branch_identifier_1;
        let depth_in_loop_1 = self.rmq.tree[node_1].depth - loop_1_length;

        debug_assert!(loop_1_length + depth_in_loop_1 == self.rmq.tree[node_1].depth);
        debug_assert!(self.rmq.tree[node_1].depth < self.rmq.tree[node_2].depth);

        let sample_ancestor  = self.la_query(         node_2, depth_in_loop_1);
        let shifted_ancestor = self.la_query(sample_ancestor, loop_1_length);
        let (rmq_start, rmq_end) = self.find_rmq_arguments(sample_ancestor, shifted_ancestor);
        let shift_extension = self.rmq.rmq.query(rmq_start, rmq_end);

        let total_subtree_extension = shift_extension + depth_in_loop_1;
        if total_subtree_extension < self.rmq.tree[node_2].depth {
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
        let k = self.k;

        debug_assert!(k < branch_identifier_1 && branch_identifier_1 < 2 * k);
        debug_assert!(k < branch_identifier_2 && branch_identifier_2 < 2 * k);

        let mut loop_1_length   = self.rmq.tree[node_1].depth + k - branch_identifier_1;
        let mut depth_in_loop_1 = self.rmq.tree[node_1].depth - loop_1_length;
        let mut loop_2_length   = self.rmq.tree[node_2].depth + k - branch_identifier_2;
        let mut depth_in_loop_2 = self.rmq.tree[node_2].depth - loop_2_length;

        let (rmq_start, rmq_end) = self.find_rmq_arguments(node_1, node_2);
        let initial_extension = self.rmq.rmq.query(rmq_start, rmq_end);

        if initial_extension < self.rmq.tree[node_1].depth && initial_extension < self.rmq.tree[node_2].depth {
            let ordering = self.compare_ancestors(node_1, node_2, initial_extension);
            return (initial_extension, ordering);
        }

        if initial_extension >= self.rmq.tree[node_1].depth && initial_extension >= self.rmq.tree[node_2].depth {
            return (k, std::cmp::Ordering::Equal);
        }

        let mut reverse = false;
        if initial_extension >= self.rmq.tree[node_2].depth {
            reverse = true;
            std::mem::swap(&mut node_1,              &mut node_2);
            std::mem::swap(&mut branch_identifier_1, &mut branch_identifier_2);
            std::mem::swap(&mut loop_1_length,       &mut loop_2_length);
            std::mem::swap(&mut depth_in_loop_1,     &mut depth_in_loop_2);
        }

        // (initial_extension >= self.rmq.tree[node_1].depth and initial_extension <  self.rmq.tree[node_2].depth)

        let sample_ancestor  = self.find_canonic_node(self.la_query(         node_2, depth_in_loop_1));
        let shifted_ancestor = self.find_canonic_node(self.la_query(sample_ancestor, loop_1_length));
        let (rmq_start, rmq_end) = self.find_rmq_arguments(sample_ancestor, shifted_ancestor);
        let shift_extension = self.rmq.rmq.query(rmq_start, rmq_end);
        
        // longer loop === shifted_ancestor
        let mut ordering = self.compare_ancestors(sample_ancestor, shifted_ancestor, shift_extension);
        if reverse {
            ordering = ordering.reverse();
        }

        (shift_extension + depth_in_loop_1, ordering)
    }

    #[inline]
    fn la_query(&self, node: usize, climb: usize) -> usize {
        self.laq.query(&self.rmq.tree, node, climb)
    }

    #[inline]
    fn find_rmq_arguments(&self, node_1: usize, node_2: usize) -> (usize, usize) {
        let name_1 = self.rmq.tree[node_1].name;
        let name_2 = self.rmq.tree[node_2].name;
        let smaller_name = name_1.min(name_2);
        let rmq_start = smaller_name + 1;
        let rmq_end = name_1 + name_2 - smaller_name;
        (rmq_start, rmq_end)
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
        let parent_1_name = self.rmq.tree[parent_1].name;
        let parent_2_name = self.rmq.tree[parent_2].name;
        parent_1_name.cmp(&parent_2_name)
    }

    #[inline]
    fn find_branch_identifier(&self, node: usize) -> usize {
        let mut sbwt_index = self.rmq.tree[node].id;
        if sbwt_index < self.offset {
            self.rmq.sbwt_1_mapping.branch_identifier[sbwt_index]
        } else {
            sbwt_index -= self.offset;
            self.rmq.sbwt_2_mapping.branch_identifier[sbwt_index]
        }
    }

    #[inline]
    fn find_canonic_node(&self, node: usize) -> usize {
        let mut sbwt_index = self.rmq.tree[node].id;
        if sbwt_index < self.offset {
            self.rmq.sbwt_1_mapping.sbwt_to_tree[sbwt_index]
        } else {
            sbwt_index -= self.offset;
            self.rmq.sbwt_2_mapping.sbwt_to_tree[sbwt_index]
        }
    }

    #[inline]
    fn find_canonic_node_and_branch_identifier(&self, node: usize) -> (usize, usize) {
        let mut sbwt_index = self.rmq.tree[node].id;
        if sbwt_index < self.offset {
            (
                self.rmq.sbwt_1_mapping.sbwt_to_tree[sbwt_index],
                self.rmq.sbwt_1_mapping.branch_identifier[sbwt_index]
            )
        } else {
            sbwt_index -= self.offset;
            (
                self.rmq.sbwt_2_mapping.sbwt_to_tree[sbwt_index],
                self.rmq.sbwt_2_mapping.branch_identifier[sbwt_index]
            )
        }
    }
}

