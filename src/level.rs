// https://codeforces.com/blog/entry/126580
// https://judge.yosupo.jp/submission/242808

use crate::dctree::Node;

#[derive(Debug)]
pub struct NodeInfo {
    pub ladder_start: usize,
    pub ladder_length: usize,
    pub euler_index: i64,
    pub deepest_leaf: usize,
}

#[derive(Debug)]
pub struct LevelAncestorSupport {
    pub jump_pointers: Vec<(usize, usize)>,
    pub ladder_values: Vec<usize>,
    pub auxiliary_info: Vec<NodeInfo>,
}

impl LevelAncestorSupport {
    pub fn new(tree: &[Node]) -> Self {
        let tree_size = tree.len();

        let mut jump_pointers  = Vec::<(usize, usize)>::with_capacity(2 * tree_size);
        let mut ladder_values  = Vec::<usize>::with_capacity(3 * tree_size);
        let mut auxiliary_info = Vec::<NodeInfo>::with_capacity(tree_size);

        let mut stack = Vec::<usize>::new();

        let mut euler_index: i64 = 0;
        #[inline(always)]
        fn push_jump_pointers(jump_pointers: &mut Vec<(usize, usize)>, stack: &[usize], euler_index: i64) {
            // let largest_power_dividing = euler_index.isolate_lowest_one();
            let largest_power_dividing = euler_index ^ (euler_index & (euler_index - 1));
            let index_in_stack_1 = 0.max(stack.len() as i64 - 1 - 2 * largest_power_dividing);
            let index_in_stack_2 = 0.max(stack.len() as i64 - 1 - 4 * largest_power_dividing);
            let jump_1 = stack[index_in_stack_1 as usize];
            let jump_2 = stack[index_in_stack_2 as usize];
            jump_pointers.push((jump_1, jump_2));
        }

        for (node_index, node) in tree.iter().enumerate() {
            let mut previous_deepest_leaf = tree_size;
            while let Some(&top) = stack.last() && top >= node.parent {
                if previous_deepest_leaf == tree_size {
                    previous_deepest_leaf = auxiliary_info[top].deepest_leaf;
                    // We have already pushed the jump pointers for this instance of the node in
                    // the euler tour.
                } else {
                    let previous_depth = tree[previous_deepest_leaf].depth;
                    let current_deepest_leaf = auxiliary_info[top].deepest_leaf;
                    let current_depth = tree[current_deepest_leaf].depth;
                    if previous_depth > current_depth {
                        auxiliary_info[top].deepest_leaf = previous_deepest_leaf;
                    }
                    previous_deepest_leaf = auxiliary_info[top].deepest_leaf;

                    push_jump_pointers(&mut jump_pointers, &stack, euler_index);
                }
                
                if top > node.parent {
                    stack.pop();
                } else {
                    break;
                }
                euler_index += 1;
            }

            debug_assert!(stack.is_empty() || *stack.last().unwrap() == node.parent);

            stack.push(node_index);
            euler_index += 1;

            auxiliary_info.push(NodeInfo {
                ladder_start: 0,
                ladder_length: 0,
                euler_index,
                deepest_leaf: node_index,
            });

            push_jump_pointers(&mut jump_pointers, &stack, euler_index);
        }

        let mut previous_deepest_leaf = tree_size;
        while let Some(&top) = stack.last() {
            if previous_deepest_leaf == tree_size {
                previous_deepest_leaf = auxiliary_info[top].deepest_leaf;
            } else {
                let previous_depth = tree[previous_deepest_leaf].depth;
                let current_deepest_leaf = auxiliary_info[top].deepest_leaf;
                let current_depth = tree[current_deepest_leaf].depth;
                if previous_depth > current_depth {
                    auxiliary_info[top].deepest_leaf = previous_deepest_leaf;
                }
                previous_deepest_leaf = auxiliary_info[top].deepest_leaf;
                push_jump_pointers(&mut jump_pointers, &stack, euler_index);
            }
            stack.pop();
            euler_index += 1;
        }
 
        for (node_index, node) in tree.iter().enumerate() {
            let parent = node.parent;
            let parent_deepest_leaf = auxiliary_info[parent].deepest_leaf;
            let deepest_leaf = auxiliary_info[node_index].deepest_leaf;
            if node_index == 0 || deepest_leaf != parent_deepest_leaf {
                // Each deepest leaf will be initialised once so there's no need to allocate a
                // buffer for each ladder.

                let node_info = &mut auxiliary_info[deepest_leaf];
                node_info.ladder_start = ladder_values.len();
                let deepest_leaf_depth = tree[deepest_leaf].depth;
                let ladder_length = ((deepest_leaf_depth - node.depth) * 2).min(deepest_leaf_depth); 
                node_info.ladder_length = ladder_length;

                let mut ladder_rung = tree[deepest_leaf].parent;
                for i in 0..ladder_length {
                    ladder_values.push(ladder_rung);
                    ladder_rung = tree[ladder_rung].parent;
                }
            }
        }

        LevelAncestorSupport {
            jump_pointers,
            ladder_values,
            auxiliary_info,
        }
    }

    pub fn query(&self, tree: &[Node], node: usize, climb: usize) -> usize {
        match climb {
            0 => node,
            1 => tree[node].parent,
            2 => tree[tree[node].parent].parent,
            _ => {
                let target_power_of_2 = 1_i64 << msb(climb / 3);
                // println!("target_power_of_2: {}, climb: {}", target_power_of_2, climb);
                let euler_index_of_node = self.auxiliary_info[node].euler_index;
                let euler_index = (euler_index_of_node & (-target_power_of_2)) | target_power_of_2;

                // Euler index is 1 based.
                let (jump_1, jump_2) = self.jump_pointers[euler_index as usize - 1];
                let leaf = if tree[node].depth - tree[jump_2].depth <= climb {
                    self.auxiliary_info[jump_2].deepest_leaf
                } else {
                    self.auxiliary_info[jump_1].deepest_leaf
                };

                let ladder_start = self.auxiliary_info[leaf].ladder_start;
                let rungs = climb + tree[leaf].depth - tree[node].depth;
                debug_assert!(rungs <= self.auxiliary_info[leaf].ladder_length);

                // ladder values do not contain the node itself
                self.ladder_values[ladder_start + rungs - 1]
            }
        }
    }
}

#[inline(always)]
pub fn msb(value: usize) -> u32 {
    usize::BITS - 1 - value.leading_zeros()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;
    use bitvec::prelude::*;
    use sbwt::BuildByLibsais;
    use sbwt::SbwtIndex;
    use sbwt::SubsetSeq;
    use sbwt::SubsetMatrix;

    fn level_ancestor_query(tree: &[Node], mut node: usize, mut climb: usize) -> usize {
        while climb > 0 {
            node = tree[node].parent;
            climb -= 1;
        }
        node
    }

    fn test_level_ancestor(tree: &[Node]) {
        let la = LevelAncestorSupport::new(tree);
        for (node_index, node) in tree.iter().enumerate() {
            let max_depth = node.depth;
            for climb in 0..max_depth {
                let invariant = level_ancestor_query(tree, node_index, climb);
                let found     = la.query            (tree, node_index, climb);
                assert_eq!(invariant, found, "node_index: {}; climb: {}", node_index, climb);
            }
        }
    }

    #[test]
    fn test_1() {
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
            b"ACACCACCCACCCC",
            b"ACACCAGGGAGGGG",
            b"ACAGGAGGGAGGGG",
        ];
        let (mut sbwt_2, _) = BitPackedKmerSortingMem::new_from_vecs(&seq_2, k)
            .build_select_support(true)
            .run();
        let sbwts = &[
            &SbwtIndexVariant::SubsetMatrix(sbwt_1),
            &SbwtIndexVariant::SubsetMatrix(sbwt_2),
        ];
        let tp = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        tp.scope(|s| {
            s.spawn(|_| {
                let decomposition = crate::dctree::TreeDecompositionRMQ::new(1, sbwts);
                test_level_ancestor(&decomposition.tree);
            });
        });
    }

    fn node(parent: usize) -> Node {
        Node { id: 0, depth: 0, parent, name: 0 }
    }

    #[test]
    fn test_2() {
        let mut tree = vec![
            node(0), // 0
            node(0), // 1
            node(1), // 2
            node(2), // 3
            node(3), // 4
            node(4), // 5
            node(3), // 6
            node(6), // 7
            node(7), // 8
            node(8), // 9
            node(9), // 10
        ];
        for i in 1..tree.len() {
            tree[i].depth = tree[tree[i].parent].depth + 1;
        }
        test_level_ancestor(&tree);
    }
}

