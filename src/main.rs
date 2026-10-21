#![allow(unused)]

use sbwt::BuildByLibsais;
use sbwt::SbwtIndex;
use sbwt::SubsetSeq;
use sbwt::SubsetMatrix;

struct Node {
    id: usize,
    depth: usize,
    parent: usize,
    name: usize,
}

struct RangeMinimumSupport;
struct LongestCommonExtensionSupport;
fn longest_common_extension() -> LongestCommonExtensionSupport {
    // let tree: [Node; ..] = ;
    //     identity         = [ (merged common) root, sbwt1 nodes in preorder, sbwt2 nodes in preorder ];
    // let loop_anchors     = [ 0, index in identity ];
    // let length           = [..]; // to identify dummy k-mers?
    LongestCommonExtensionSupport
}

// #[inline(always)]
// fn parent2(tree: &[Node], node: usize) -> usize {
//     tree[tree[node].parent].parent
// }

#[inline]
fn climb(tree: &[Node], node: usize, offset: usize) -> usize {
    match offset {
        0 => node,
        1 => tree[node].parent,
        _ => tree[tree[node].parent].parent,
    }
}

/// Expects the names to be in [1..=alpha].
fn radix_pass(
    tree: &[Node],
    source_buffer: &[usize],
    target_buffer: &mut [usize],
    offset: usize,
    bucket_buffer: &mut Vec<usize>,
    alpha: usize
) {
    bucket_buffer.clear();
    bucket_buffer.reserve(alpha + 1);
    bucket_buffer.extend(std::iter::repeat_n(0, alpha + 1));

    for &node in source_buffer {
        let offset_node = climb(tree, node, offset);
        let name = tree[offset_node].name;
        bucket_buffer[name] += 1;
    }

    for i in 1..alpha {
        bucket_buffer[i] += bucket_buffer[i - 1]
    }

    for &node in source_buffer {
        let offset_node = climb(tree, node, offset);
        let name = tree[offset_node].name - 1;
        let target_index = bucket_buffer[name];
        bucket_buffer[name] += 1;
        target_buffer[target_index] = node;
    }
}

const BUFFER_SPACE: usize = 1;

// alpha - alphabet size
#[allow(clippy::ptr_arg)]
#[allow(clippy::only_used_in_recursion)]
fn sort_tree(
    tree: &mut [Node],
    order: &mut [usize],
    bucket_buffer: &mut Vec<usize>,
    alpha: usize,
    depth: usize,
    k: usize
) {
    let n = tree.len();

    // Step 0: count the nodes at the different depths.
    let mut modulo_counts = [0_usize; 3];
    tree[0].depth = 0;
    for i in 1..n {
        tree[i].depth = tree[tree[i].parent].depth + 1;
        let modulo = tree[i].depth % 3;
        modulo_counts[modulo] += 1;
    }

    // Choose which nodes not to include in the sample such that the number of sampled nodes is
    // guaranteed to decrease by at least 2/3.
    let target_count = n / 3;
    let mut residual_modulo = 0;
    #[allow(clippy::needless_range_loop)]
    for modulo in 0..3 {
        if modulo_counts[modulo] > target_count {
            residual_modulo = modulo;
        }
    }

    let mut residual_count = modulo_counts[residual_modulo];
    let mut sample_count = n - residual_count;

    // Include the root node in the sample either way.
    if residual_modulo == 0 {
        sample_count += 1;
        residual_count -= 1;
    }

    // Step 0.1: construct the sample.
    let mut tree_to_sample_tree_mapping = vec![0_usize; n];
    let mut sample_tree         = Vec::<Node>::with_capacity(sample_count);
    let mut sample_indices      = Vec::<usize>::with_capacity(sample_count);
    // This buffer will be reused later; reserve some more space.
    let mut sample_suffix_array = vec![0_usize; sample_count + BUFFER_SPACE];
    sample_suffix_array.pop();
    let mut residual_indices    = Vec::<usize>::with_capacity(residual_count);

    sample_indices.push(0);
    sample_tree.push(Node {
        id: 0,
        depth: 0,
        parent: 0,
        name: 0,
    });

    for i in 1..n {
        if tree[i].depth % 3 != residual_modulo {
            let parent = climb(tree, i, 2);
            let parent_in_sample = tree_to_sample_tree_mapping[parent];
            tree_to_sample_tree_mapping[i] = sample_tree.len();
            sample_tree.push(Node {
                id: i,
                depth: 0,
                parent: parent_in_sample,
                name: 0,
            });
            sample_indices.push(i);
        } else {
            residual_indices.push(i);
        }
    }

    assert_eq!(sample_indices.len(), sample_count);
    assert_eq!(sample_tree.len(), sample_count);
    assert_eq!(residual_indices.len(), residual_count);

    // Step 1: sort the samples.
    radix_pass(tree, &sample_indices,      &mut sample_suffix_array, 2, bucket_buffer, alpha);
    radix_pass(tree, &sample_suffix_array, &mut sample_indices,      1, bucket_buffer, alpha);
    radix_pass(tree, &sample_indices,      &mut sample_suffix_array, 0, bucket_buffer, alpha);

    // Step 1.1: generate new names for the samples
    let null_alpha = alpha + 1;
    let mut tuple_name: (usize, usize, usize) = (null_alpha, null_alpha, null_alpha);
    let mut name = 0;
    for i in 0..sample_count {
        let index_in_tree = sample_suffix_array[i];
        let sampled_node_tuple_name = (
            climb(tree, index_in_tree, 2),
            climb(tree, index_in_tree, 1),
            climb(tree, index_in_tree, 0),
        );
        if tuple_name != sampled_node_tuple_name {
            name += 1;
            tuple_name = sampled_node_tuple_name;
        }
        let index_in_sample_tree = tree_to_sample_tree_mapping[index_in_tree];
        sample_tree[index_in_sample_tree].name = name;
        sample_indices[i] = index_in_sample_tree;
    }

    // Step 1.2: check whether the names are unique.
    if name < sample_tree.len() {
        // TODO: check whether depth is enough to sort the first k characters.
        let alpha = name;
        sort_tree(&mut sample_tree, &mut sample_indices, bucket_buffer, alpha, depth * 3, k);
        #[allow(clippy::needless_range_loop)]
        for i in 0..sample_count {
            let index_in_sample_tree = sample_indices[i];
            sample_tree[index_in_sample_tree].name = i;
        }
    }

    // Step 2: sort residual indices.
    //
    // Radix sort using the new names from the sample tree.
    // note(mk): (?) reuse one of the fields of Node to skip using the mapping a second time.
    let mut residual_buckets = sample_suffix_array; 
    residual_buckets.iter_mut().for_each(|value| *value = 0);
    residual_buckets.push(0);
    for &residual_index in &residual_indices {
        let parent_in_tree   = tree[residual_index].parent;
        let parent_in_sample = tree_to_sample_tree_mapping[parent_in_tree];
        let parent_name = sample_tree[parent_in_sample].name;
        residual_buckets[parent_name] += 1;
    }
    for i in 1..=sample_count {
        residual_buckets[i] += residual_buckets[i - 1];
    }
    for &residual_index in &residual_indices {
        let parent_in_tree   = tree[residual_index].parent;
        let parent_in_sample = tree_to_sample_tree_mapping[parent_in_tree];
        let parent_name = sample_tree[parent_in_sample].name - 1;
        let target_index = residual_buckets[parent_name];
        residual_buckets[parent_name] += 1;
        order[target_index] = residual_index;
    }
    radix_pass(tree, &order[..residual_count], &mut residual_indices, 0, bucket_buffer, alpha);

    // Step 3: merge.
    let mut residual_i = 0;
    let mut sample_i = 0;
    while residual_i < residual_count && sample_i < sample_count {
        residual_i += 1;
        sample_i += 1;
    }

    // TODO: ...
}

struct LevelAncestorSupport;

// (Idea) construct LCP out of [order] and [tree]:
// - construct "Phi" array - i.e. for each node which other node it should be compared to by
//   scanning [order]
// - traverse the [tree] array back to front: we want the children first and then the parents;
//   when the LCP value for a child has been calculated set the "potential" LCP value of the parent
//   to be the maximum of its current "potential" LCP value (initially 0) and LCP[child]-1 
// - when calculating the LCP value for the current node use a level ancestor query to jump
//   at the parent node at which the potential LCP ends; do the same with the node from Phi;
//   find the final LCP value by traversing the tree upwards.

fn main() {
    let seq_1 = vec![b"TCGCAT".to_vec(), b"ACTCAA".to_vec()];
    let (sbwt_1, _) = BuildByLibsais::new_from_vecs(&seq_1, 3)
        .build_select_support(true)
        .run();
    let seq_2 = vec![b"GCTGGA".to_vec(), b"ACTCAA".to_vec()];
    let (sbwt_2, _) = BuildByLibsais::new_from_vecs(&seq_2, 3)
        .build_select_support(true)
        .run();
    // TODO: ...
}

