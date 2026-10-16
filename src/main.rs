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

fn sort_sbwt_graphs() {
    // let identity         = [ (merged common) root, sbwt1 nodes in preorder, sbwt2 nodes in preorder ];
    // let loop_anchors     = [ 0, index in identity ];
    // let tree: [Node; ..] = ;
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
        bucket_buffer[name + 1] += 1;
    }

    for i in 1..alpha {
        bucket_buffer[i] += bucket_buffer[i - 1]
    }

    for &node in source_buffer {
        let offset_node = climb(tree, node, offset);
        let name = tree[offset_node].name;
        let target_index = bucket_buffer[name];
        bucket_buffer[name] += 1;
        target_buffer[target_index] = node;
    }
}

//
// alpha - alphabet size
//
fn sort_tree(
    tree: &mut [Node],
    order: &mut [usize],
    bucket_buffer: &mut Vec<usize>,
    alpha: usize,
    depth: usize,
    k: usize
) {

    // 
    // // n := tree.len();
    //
    // // step 0: construct the sample (calculate depth?) count modulos.
    // let tree_to_sampled_set_mapping : Vec<usize> = [n];
    // //  contains mapping to the original tree.
    // //  vvvvvvvvvvvvv 
    // let sampled_nodes               : Vec<Node> = vec![];
    // let sampled_set                 : Vec<usize> = [2/3n];
    // let sampled_set_swap            : Vec<usize> = [2/3n];
    //
    // let residual_set                : Vec<usize> = [1/3n];
    //
    // // - [residual_set] is a subset of [0..n)
    // // - for i in residual_set we can find the parent j of i in [tree]
    // // - from [tree_to_sampled_set_mapping] and j we can find the index j' in [sampled_set]
    //
    // // step 1: sort the samples
    // radix_pass(tree, sampled_set, sampled_set_swap, offset: 2);
    // radix_pass(tree, sampled_set_swap, sampled_set, offset: 1);
    // radix_pass(tree, sampled_set, sampled_set_swap, offset: 0);
    //
    // // - [sampled_set] contains (almost) sorted indices in [tree]
    //
    // // step 1.1: generate new names for the samples
    // let mut name = 0;
    //
    // // - in the order induced by [sampled_set] iterate through [sampled_nodes]
    // // ! pass through the mapping to get the correct index
    // //   - if the 3-tuple of characters of this node is different from the previous:
    // //     1. assign it the current name (the first 3-tuple is null)
    // //     2. increment the name
    // //   - set the name of this node
    // 
    //
    // // TODO: rework this:
    //
    // // step 1.2: check whether the names are unique and generate the order of the nodes
    // if name < sampled_nodes.len() {
    //     // the names are not unique
    //     // alpha = name
    //     let order_of_sampled_nodes = order_nodes(
    //         &mut sampled_nodes,
    //         &mut sampled_node_indices,
    //         name,
    //         depth + 1,
    //         k
    //     );
    //     // store the unique names for the sampled_nodes
    // } else {
    //     // all names are unique
    //     // [sampled_node_indices_swap] contains indices into [tree]
    //     // we want [sampled_node_indices] to contain 
    // }
    //
    // // step 2: sort [residual_set]
    // // move [sampled_set_swap] to [residual_set_buckets]
    // radix_pass(tree, residual_set, residual_set_buckets, offset: 0);
    //
    // 
    // TODO: finish sketching out the algorithm
    //
}

fn main() {

}

