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

#[inline(always)]
fn parent2(tree: &[Node], node: usize) -> usize {
    tree[tree[node].parent].parent
}

#[inline(always)]
fn parent3(tree: &[Node], node: usize) -> usize {
    tree[parent2(tree, node)].parent
}

#[inline]
fn climb(tree: &[Node], node: usize, offset: usize) -> usize {
    match offset {
        0 => node,
        1 => tree[node].parent,
        2 => parent2(tree, node),
        _ => parent3(tree, node),
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
    bucket_buffer.clear();
    bucket_buffer.reserve(alpha);

    // 
    // // n := tree.len();
    //
    // // step 0: construct the sample (calculate depth?)
    // // count modulos.
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
    // radix_pass(tree, sampled_set_swap, sampled_set, parent_offset: 2);
    // radix_pass(tree, sampled_set, sampled_set_swap, parent_offset: 1);
    // radix_pass(tree, sampled_set_swap, sampled_set, parent_offset: 0);
    // // copy [sampled_set] to [sampled_set_swap]
    //
    // // - [sampled_set] contains (almost) sorted indices in [tree]
    //
    // // step 1.1: generate new names for the samples
    // let mut name = 0;
    //
    // // - in the order induced by [sampled_set] iterate through [sampled_nodes]
    // //   - if the 3-tuple of characters of this node is different from the previous:
    // //     1. assign it the current name (the first 3-tuple is null)
    // //     2. increment the name
    // //   - set the name
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
    // 
    // TODO: finish sketching out the algorithm
    //
}

fn main() {

}

