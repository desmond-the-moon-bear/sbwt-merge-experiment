#![allow(unused)]

use sbwt::BuildByLibsais;
use sbwt::SbwtIndex;
use sbwt::SubsetSeq;
use sbwt::SubsetMatrix;

#[derive(Default, Debug)]
struct Node {
    id: usize,
    depth: usize,
    parent: usize,
    name: usize,
}

struct RangeMinimumSupport;
struct LongestCommonExtensionSupport;
fn longest_common_extension<SS1, SS2>(sbwt_1: &SbwtIndex<SS1>, sbwt_2: &SbwtIndex<SS2>) -> LongestCommonExtensionSupport
where
    SS1: SubsetSeq,
    SS2: SubsetSeq,
{
    // let tree: [Node; ..] = ;
    //     identity         = [ (merged common) root, sbwt1 nodes in preorder, sbwt2 nodes in preorder ];
    // let loop_anchors     = [ 0, index in identity ];
    // let length           = [..]; // to identify dummy k-mers?

    // Merge the root of the two SBWT indices.
    let sbwt_1_size = sbwt_1.n_sets();
    let sbwt_2_size = sbwt_2.n_sets();
    let tree_size = sbwt_1_size + sbwt_2_size - 1;
    let mut tree = Vec::<Node>::with_capacity(tree_size);
    tree.push(Node::default());
    tree[0].name = 1;

    let mut stack = Vec::<(usize, usize, usize)>::new();
    enumerate_tree_from_sbwt(sbwt_1, &mut stack, &mut tree, 0);
    enumerate_tree_from_sbwt(sbwt_2, &mut stack, &mut tree, sbwt_1_size);

    let mut order = vec![0_usize; tree.len()];
    let mut bucket_buffer = vec![];
    sort_tree(&mut tree, &mut order, &mut bucket_buffer, 5, 1, sbwt_1.k());

    let mut kmer_buffer = vec![0_u8; sbwt_1.k()];
    for &tree_index in &order {
        let node_index = tree[tree_index].id;
        kmer_buffer.clear();
        if node_index < sbwt_1_size {
            sbwt_1.push_kmer_to_vec(node_index, &mut kmer_buffer);
            let kmer = str::from_utf8(&kmer_buffer).unwrap();
            println!("[1:{:05}] {}", node_index, kmer);
        } else {
            let node_index = node_index - sbwt_1_size;
            sbwt_2.push_kmer_to_vec(node_index, &mut kmer_buffer);
            let kmer = str::from_utf8(&kmer_buffer).unwrap();
            println!("[2:{:05}] {}", node_index, kmer);
        }
    }

    LongestCommonExtensionSupport
}

fn enumerate_tree_from_sbwt<SS: SubsetSeq>(
    sbwt: &SbwtIndex<SS>,
    stack: &mut Vec<(usize, usize, usize)>,
    tree: &mut Vec<Node>,
    offset: usize
) {
    let alphabet = sbwt.alphabet().iter();
    stack.push((0, 0, 1));
    while let Some((sbwt_index, tree_parent, name)) = stack.pop() {
        let tree_index = if sbwt_index != 0 {
            let index = tree.len();
            tree.push(Node {
                id: sbwt_index + offset,
                depth: 0,
                parent: tree_parent,
                name,
            });
            index
        } else {
            0
        };

        for (character_index, &character) in alphabet.clone().enumerate().rev() {
            let option_range = sbwt.search_from(sbwt_index..sbwt_index + 1, &[character]);
            if let Some(range) = option_range {
                stack.push((range.start, tree_index, character_index + 1));
            }
        }
    }
}

#[inline]
fn climb(tree: &[Node], node: usize, offset: usize) -> usize {
    match offset {
        0 => node,
        1 => tree[node].parent,
        2 => tree[tree[node].parent].parent,
        _ => tree[tree[tree[node].parent].parent].parent,
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

const EXTRA_BUFFER_SPACE: usize = 1;

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
    modulo_counts[0] += 1;
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
    let mut sample_suffix_array = vec![0_usize; sample_count + EXTRA_BUFFER_SPACE];
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
            let parent = climb(tree, i, 3);
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
            tree[climb(tree, index_in_tree, 0)].name,
            tree[climb(tree, index_in_tree, 1)].name,
            tree[climb(tree, index_in_tree, 2)].name,
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
    if name < sample_tree.len() && depth < k {
        let alpha = name;
        sort_tree(&mut sample_tree, &mut sample_indices, bucket_buffer, alpha, depth * 3, k);
        #[allow(clippy::needless_range_loop)]
        for i in 0..sample_count {
            let index_in_sample_tree = sample_indices[i];
            sample_tree[index_in_sample_tree].name = i + 1;
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
    assert_eq!(residual_count + sample_count, order.len());
    let mut residual_i = 0;
    let mut sample_i = 0;
    let mut write_i = 0;
    while residual_i < residual_count && sample_i < sample_count {
        let current_residual_node = residual_indices[residual_i];
        let index_in_sample_tree  = sample_indices[sample_i];
        let current_sample_node   = sample_tree[index_in_sample_tree].id;
        let current_sample_modulo = tree[current_sample_node].depth % 3;

        if (current_sample_modulo + 2) % 3 == residual_modulo {
            let sample_parent = climb(tree, current_residual_node, 2);
            let sample_parent_in_sample_tree = tree_to_sample_tree_mapping[sample_parent];
            let sample_parent_name = sample_tree[sample_parent_in_sample_tree].name;
            let residual_node_name = (
                tree[current_residual_node].name,
                tree[climb(tree, current_residual_node, 1)].name,
                sample_parent_name,
                current_residual_node,
            );

            let sample_parent = climb(tree, current_sample_node, 2);
            let sample_parent_in_sample_tree = tree_to_sample_tree_mapping[sample_parent];
            let sample_parent_name = sample_tree[sample_parent_in_sample_tree].name;
            let sample_node_name = (
                tree[current_sample_node].name,
                tree[climb(tree, current_sample_node, 1)].name,
                sample_parent_name,
                current_sample_node,
            );

            if residual_node_name < sample_node_name {
                order[write_i] = current_residual_node;
                residual_i += 1;
            } else {
                order[write_i] = current_sample_node;
                sample_i += 1;
            }
        } else {
            let sample_parent = climb(tree, current_residual_node, 1);
            let sample_parent_in_sample_tree = tree_to_sample_tree_mapping[sample_parent];
            let sample_parent_name = sample_tree[sample_parent_in_sample_tree].name;
            let residual_node_name = (
                tree[current_residual_node].name,
                sample_parent_name,
                current_residual_node,
            );

            let sample_parent = climb(tree, current_sample_node, 1);
            let sample_parent_in_sample_tree = tree_to_sample_tree_mapping[sample_parent];
            let sample_parent_name = sample_tree[sample_parent_in_sample_tree].name;
            let sample_node_name = (
                tree[current_sample_node].name,
                sample_parent_name,
                current_sample_node,
            );

            if residual_node_name < sample_node_name {
                order[write_i] = current_residual_node;
                residual_i += 1;
            } else {
                order[write_i] = current_sample_node;
                sample_i += 1;
            }
        }

        write_i += 1;
    }
    while residual_i < residual_count {
        let current_residual_node = residual_indices[residual_i];
        order[write_i] = current_residual_node;
        residual_i += 1;
        write_i += 1;
    }
    while sample_i < sample_count {
        let index_in_sample_tree  = sample_indices[sample_i];
        let current_sample_node   = sample_tree[index_in_sample_tree].id;
        order[write_i] = current_sample_node;
        sample_i += 1;
        write_i += 1;
    }
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
    let seq_1 = vec![
        b"TCGCCCTCGCAT".to_vec(),
        b"ACTCTGACTCAA".to_vec()
    ];
    let (sbwt_1, _) = BuildByLibsais::new_from_vecs(&seq_1, 10)
        .build_select_support(true)
        .run();

    let seq_2 = vec![
        b"GCTGGTGCTTGA".to_vec(),
        b"ACTCTGACTCAA".to_vec()
    ];
    let (sbwt_2, _) = BuildByLibsais::new_from_vecs(&seq_2, 10)
        .build_select_support(true)
        .run();

    let lce = longest_common_extension(&sbwt_1, &sbwt_2);
}

