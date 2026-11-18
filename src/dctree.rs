use crate::rmq::RangeMinimumQuery as RMQ;
use crate::level::LevelAncestorSupport as LAQ;
use crate::lce::LongestCommonExtension as LCE;

use bitvec::prelude::*;
use sbwt::BuildByLibsais;
use sbwt::BitPackedKmerSortingMem;
use sbwt::SubsetSeq;
use sbwt::SubsetMatrix;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

#[derive(Default, Debug)]
pub struct Node {
    pub id: usize,
    pub depth: usize,
    pub parent: usize,
    pub name: usize,
}

#[derive(Default, Debug)]
pub struct MappingInformation {
    pub max_sbwt_size: usize,
    pub sbwt_offsets: Vec<usize>,
    pub sbwt_to_tree: Vec<usize>,
    // (type:1) [0..=k]      - branch originally connected to the root
    //                         or of a loop with length greater than k
    // (type:2) [k+1..=2k-1] - node in an unrolled loop
    // (type:3) [2k..]       - in the subtree of a node in an unrolled loop;
    //                         points to the root of the subtree
    pub branch_identifier: Vec<usize>,
}

impl MappingInformation {
    fn new(sbwts: &[&SbwtIndexVariant]) -> Self {
        let mut sbwt_offsets = vec![0; sbwts.len() + 1];
        let mut max_sbwt_size = 0;
        let mut current_size = 0;
        let mut total_size = 0;

        for index in 0..sbwts.len() {
            current_size  = sbwts[index].n_sets();
            max_sbwt_size = max_sbwt_size.max(current_size);
            sbwt_offsets[index + 1] = sbwt_offsets[index] + current_size;
            total_size   += current_size;
        }

        let sbwt_to_tree = vec![0_usize; total_size];
        let branch_identifier = vec![0_usize; total_size];

        Self {
            max_sbwt_size,
            sbwt_offsets,
            sbwt_to_tree,
            branch_identifier,
        }
    }

    pub fn set(&mut self, colex_rank: usize, origin: usize, node: usize, branch_identifier: usize) {
        let mapping_index = self.sbwt_offsets[origin] + colex_rank;
        self.sbwt_to_tree[mapping_index] = node;
        self.branch_identifier[mapping_index] = branch_identifier;

    }

    pub fn node(&self, colex_rank: usize, origin: usize) -> usize {
        self.sbwt_to_tree[self.sbwt_offsets[origin] + colex_rank]
    }

    pub fn branch_id(&self, colex_rank: usize, origin: usize) -> usize {
        self.branch_identifier[self.sbwt_offsets[origin] + colex_rank]
    }

    pub fn node_and_branch_id(&self, colex_rank: usize, origin: usize) -> (usize, usize) {
        let mapping_index = self.sbwt_offsets[origin] + colex_rank;
        (self.sbwt_to_tree[mapping_index], self.branch_identifier[mapping_index])
    }

    pub fn id_to_colex_rank(&self, id: usize) -> (usize, usize) {
        let origin = id / self.max_sbwt_size;
        (id - origin * self.max_sbwt_size, origin)
    }

    pub fn id_to_canonic_node(&self, id: usize) -> usize {
        let origin = id / self.max_sbwt_size;
        let colex_rank = id - origin * self.max_sbwt_size;
        self.node(colex_rank, origin)
    }
}

pub struct TreeDecompositionRMQ {
    pub mapping: MappingInformation,
    // The names in the tree should be the inverse permutation of the order vector. That is from
    // the (sbwt -> tree) mapping we can get the range for the range minimum queries.
    pub tree: Vec<Node>,
    pub rmq: RMQ,
}

pub fn tree_decomposition_rmq(sbwts: &[&SbwtIndexVariant]) -> TreeDecompositionRMQ {
    assert!(sbwts.len() > 1);
    let k = sbwts[0].k();
    #[allow(clippy::needless_range_loop)]
    for i in 1..sbwts.len() {
        assert!(k == sbwts[i].k());
    }

    // Merge the root of the two SBWT indices.
    let mut mapping = MappingInformation::new(sbwts);
    let initial_tree_capacity = 2 * mapping.sbwt_to_tree.len();
    let mut tree = Vec::<Node>::with_capacity(initial_tree_capacity);
    tree.push(Node::default());
    tree[0].name = 1;

    let mut stack = Vec::<(usize, usize, usize)>::new();
    log::info!("emitting nodes:");
    for (index, sbwt) in sbwts.iter().enumerate() {
        let mut sbwt_size = sbwt.n_sets();
        let mut mapped = bitvec![0; sbwt_size];
        let offset = index * mapping.max_sbwt_size;
        emit_tree_from_sbwt      (sbwt, index, &mut stack, &mut tree, &mut mapping, &mut mapped, offset);
        emit_components_from_sbwt(sbwt, index, &mut stack, &mut tree, &mut mapping, &mut mapped, offset);
    }

    let tree_size = tree.len();

    let mut order = vec![0_usize; tree_size];

    let mut alphabet_lcp          = Vec::with_capacity(tree_size);
    let mut alphabet_name_lengths = Vec::with_capacity(tree_size);
    alphabet_lcp         .extend_from_slice(&[0; 5]);
    alphabet_name_lengths.push(0);
    alphabet_name_lengths.extend_from_slice(&[1; 4]);
    let mut alphabet_rmq = RMQ::new(alphabet_name_lengths, alphabet_lcp);

    log::info!("sorting the tree:");
    let mut bucket_buffer = vec![];
    let rmq = sort_tree(&mut tree, &mut order, &mut alphabet_rmq, &mut bucket_buffer, 5, 1, k);

    for i in 0..tree_size {
        tree[order[i]].name = i;
    }

    TreeDecompositionRMQ {
        mapping,
        tree,
        rmq,
    }
}

fn inedge_name(sbwt: &SbwtIndexVariant, colex_rank: usize) -> usize {
    let inlabel = sbwt.inlabel(colex_rank);
    if let Some(inlabel) = inlabel {
        match inlabel {
            b'A' => 2,
            b'C' => 3,
            b'G' => 4,
            _    => 5,
        }
    } else {
        1
    }
}

fn emit_tree_from_sbwt(
    sbwt: &SbwtIndexVariant,
    origin: usize,
    stack: &mut Vec<(usize, usize, usize)>,
    tree: &mut Vec<Node>,
    mapping: &mut MappingInformation,
    mapped: &mut BitVec,
    offset: usize
) {
    let alphabet = sbwt.alphabet();
    let k = sbwt.k();

    mapped.set(0, true);
    stack.clear();
    stack.push((0, 0, 0));
    while let Some((colex_rank, tree_parent, length)) = stack.pop() {
        let mut next_edge = tree_parent;
        if !mapped[colex_rank] {
            let index = tree.len();
            next_edge = 0;
            tree.push(Node {
                id: colex_rank + offset,
                depth: 0,
                parent: tree_parent,
                name: inedge_name(sbwt, colex_rank),
            });
            mapped.set(colex_rank, true);
            mapping.set(colex_rank, origin, index, length);
        }
        let tree_index = mapping.node(colex_rank, origin);
        while next_edge < alphabet.len() {
            let character = &[alphabet[next_edge]];
            let option_range = sbwt.search_from(colex_rank..colex_rank + 1, character);
            next_edge += 1;
            if let Some(range) = option_range {
                if next_edge < alphabet.len() {
                    stack.push((colex_rank, next_edge, length));
                }
                stack.push((range.start, tree_index, (length + 1).min(k))); 
                break;
            }
        }
    }
}

fn emit_components_from_sbwt(
    sbwt: &SbwtIndexVariant,
    origin: usize,
    stack: &mut Vec<(usize, usize, usize)>,
    tree: &mut Vec<Node>,
    mapping: &mut MappingInformation,
    mapped: &mut BitVec,
    offset: usize,
) {
    let alphabet = sbwt.alphabet();
    let k = sbwt.k();
    let sbwt_size = sbwt.n_sets();

    let mut explored = bitvec![0; sbwt_size];

    for current_colex_rank in 1..sbwt_size {
        if mapped[current_colex_rank] {
            continue;
        }

        let mut next_child = 0;
        let mut found_loop = false;
        stack.clear();

        stack.push((current_colex_rank, 0, 0));

        while let Some((colex_rank, mut next_edge, _)) = stack.pop() {
            if next_edge >= alphabet.len() {
                continue;
            }

            next_child = 0;
            while next_edge < alphabet.len() {
                let option_range = sbwt.search_from(colex_rank..colex_rank + 1, &[alphabet[next_edge]]);
                next_edge += 1;
                if let Some(range) = option_range && !explored[range.start] {
                    next_child = range.start;
                    break;
                }
            }

            if next_child > 0 {
                stack.push((colex_rank, next_edge, 0));
                if next_child == current_colex_rank {
                    found_loop = true;
                    break;
                }
                explored.set(next_child, true);
                stack.push((next_child, 0, 0));
            }
        }

        explored.set(current_colex_rank, true);

        if !found_loop {
            assert!(stack.is_empty());
            continue;
        }

        assert!(!stack.is_empty());

        // First instance of the necklace:
        let mut parent = 0;
        let loop_length = stack.len();
        let to_skip = if loop_length >= k {
            stack.len() - k + 1
        } else {
            0
        };
        for &(colex_rank, _, _) in stack.iter().skip(to_skip) {
            let tree_index = tree.len();
            tree.push(Node {
                id: colex_rank + offset,
                depth: 0,
                parent,
                name: inedge_name(sbwt, colex_rank),
            });
            parent = tree_index;
        }

        if loop_length >= k {
            emit_long_loop_branch(
                sbwt, origin,
                stack,
                tree,
                mapping,
                mapped,
                offset,
                current_colex_rank,
                parent,
            );
        } else {
            emit_short_loop_branches(
                sbwt, origin,
                stack,
                tree,
                mapping,
                mapped,
                offset,
                current_colex_rank,
                parent,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn emit_long_loop_branch(
    sbwt: &SbwtIndexVariant, origin: usize,
    stack: &mut Vec<(usize, usize, usize)>,
    tree: &mut Vec<Node>,
    mapping: &mut MappingInformation,
    mapped: &mut BitVec,
    offset: usize,
    current_sbwt_index: usize,
    mut necklace_parent: usize,
) {
    let k = sbwt.k();
    let alphabet = sbwt.alphabet();

    // Explore the subtrees of the necklace:
    stack.clear();
    stack.push((current_sbwt_index, necklace_parent, 0));
    while let Some((colex_rank, tree_parent, _)) = stack.pop() {
        let mut tree_index;
        let mut next_edge = tree_parent;

        if mapped[colex_rank] {
            tree_index = mapping.node(colex_rank, origin);
        } else {
            tree_index = tree.len();
            next_edge = 0;
            tree.push(Node {
                id: colex_rank + offset,
                depth: 0,
                parent: tree_parent,
                name: inedge_name(sbwt, colex_rank),
            });
            mapped.set(colex_rank, true);
            mapping.set(colex_rank, origin, tree_index, k);
        };

        while next_edge < alphabet.len() {
            let character = &[alphabet[next_edge]];
            let option_range = sbwt.search_from(colex_rank..colex_rank + 1, character);
            next_edge += 1;
            if let Some(range) = option_range && range.start != current_sbwt_index {
                if next_edge < alphabet.len() {
                    stack.push((colex_rank, next_edge, 0));
                }
                stack.push((range.start, tree_index, 0)); 
                break;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn emit_short_loop_branches(
    sbwt: &SbwtIndexVariant, origin: usize,
    stack: &mut Vec<(usize, usize, usize)>,
    tree: &mut Vec<Node>,
    mapping: &mut MappingInformation,
    mapped: &mut BitVec,
    offset: usize,
    current_sbwt_index: usize,
    mut necklace_parent: usize,
) {
    let k = sbwt.k();
    let alphabet = sbwt.alphabet();
    let loop_length = stack.len();

    let loop_nodes_start = tree.len();
    let loop_nodes_end = loop_nodes_start + loop_length;

    // Second instance of the necklace.
    let mut loop_branch_identifier = 0;
    for &(colex_rank, _, _) in stack.iter() {
        mapped.set(colex_rank, true);
        let tree_index = tree.len();
        if colex_rank == 212 && offset > 0 {
            println!("short loop unravel");
        }
        tree.push(Node {
            id: colex_rank + offset,
            depth: 0,
            parent: necklace_parent,
            name: inedge_name(sbwt, colex_rank),
        });
        loop_branch_identifier += 1;
        mapping.set(colex_rank, origin, tree_index, k + loop_branch_identifier);
        necklace_parent = tree_index;
    }

    // Explore the subtrees of the necklace:
    #[allow(clippy::needless_range_loop)]
    for loop_node_index in loop_nodes_start..loop_nodes_end {
        let current_sbwt_index = tree[loop_node_index].id - offset;

        stack.clear();
        stack.push((current_sbwt_index, 0, 0));
        while let Some((colex_rank, tree_parent, mut branch_identifier)) = stack.pop() {
            let mut tree_index;
            let mut next_edge = tree_parent;

            debug_assert!(colex_rank == current_sbwt_index || branch_identifier != 0);

            if mapped[colex_rank] {
                tree_index = mapping.node(colex_rank, origin);
                if colex_rank == current_sbwt_index {
                    tree_index = 0;
                    branch_identifier = 2 * k + loop_node_index;
                }
            } else {
                tree_index = tree.len();
                next_edge = 0;
                tree.push(Node {
                    id: colex_rank + offset,
                    depth: 0,
                    parent: tree_parent,
                    name: inedge_name(sbwt, colex_rank),
                });
                mapped.set(colex_rank, true);
                mapping.set(colex_rank, origin, tree_index, branch_identifier);
            };

            while next_edge < alphabet.len() {
                let character = &[alphabet[next_edge]];
                let option_range = sbwt.search_from(colex_rank..colex_rank + 1, character);
                next_edge += 1;
                if let Some(range) = option_range && !mapped[range.start] {
                    if next_edge < alphabet.len() {
                        stack.push((colex_rank, next_edge, branch_identifier));
                    }
                    debug_assert!(branch_identifier >= 2 * k);
                    stack.push((range.start, tree_index, branch_identifier)); 
                    break;
                }
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
#[allow(clippy::too_many_arguments)]
fn sort_tree(
    tree: &mut [Node],
    order: &mut [usize],
    input_names: &mut RMQ,
    bucket_buffer: &mut Vec<usize>,
    alpha: usize,
    rank: usize,
    k: usize
) -> RMQ {
    let n = tree.len();

    // The caller should allocate the entire buffer as it probably doesn't know how many names the
    // samples will divide into.
    let input_names_buffer             = input_names.values().to_vec();
    let input_name_lengths_buffer      = input_names.name_lengths().to_vec();
    let mut output_lcp_values_buffer   = input_names.swap_values(input_names_buffer);
    let mut output_name_lengths_buffer = input_names.swap_name_lengths(input_name_lengths_buffer);
    output_lcp_values_buffer.clear();
    output_name_lengths_buffer.clear();
    assert!(output_lcp_values_buffer.capacity() >= n);
    assert!(output_name_lengths_buffer.capacity() >= n);

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

    let mut sample_lcp        = Vec::<usize>::with_capacity(sample_count);
    let mut sample_name_lengths = Vec::<usize>::with_capacity(sample_count);

    // Step 1.1: generate new names for the samples
    let null_alpha = alpha + 1;
    let mut tuple_name: [usize; 3] = [null_alpha, null_alpha, null_alpha];
    let mut name = 0;
    for i in 0..sample_count {
        let index_in_tree = sample_suffix_array[i];
        let sampled_node_tuple_name = [
            tree[climb(tree, index_in_tree, 0)].name,
            tree[climb(tree, index_in_tree, 1)].name,
            tree[climb(tree, index_in_tree, 2)].name,
        ];

        if tuple_name != sampled_node_tuple_name {
            if tuple_name[0] == null_alpha {
                sample_lcp       .push(0);
                sample_name_lengths.push(0);
            } else {
                let mut lcp = 0;
                #[allow(clippy::needless_range_loop)]
                for index in 0..tuple_name.len() {
                    let previous_input_name = tuple_name[index] - 1;
                    let sampled_input_name  = sampled_node_tuple_name[index] - 1;

                    if previous_input_name == sampled_input_name {
                        let length = input_names.name_lengths()[sampled_input_name];
                        lcp += length;
                        if length < rank {
                            break;
                        }
                    } else {
                        lcp += input_names.query(previous_input_name + 1, sampled_input_name);
                        break;
                    }
                }
                let mut name_length = 0;
                #[allow(clippy::needless_range_loop)]
                for index in 0..tuple_name.len() {
                    let sampled_input_name = sampled_node_tuple_name[index] - 1;
                    let length = input_names.name_lengths()[sampled_input_name];
                    name_length += length;
                    if length < rank {
                        break;
                    }
                }
                sample_lcp         .push(lcp);
                sample_name_lengths.push(name_length);
            }

            name += 1;
            tuple_name = sampled_node_tuple_name;
        }
        let index_in_sample_tree = tree_to_sample_tree_mapping[index_in_tree];
        sample_tree[index_in_sample_tree].name = name;
        sample_indices[i] = index_in_sample_tree;
    }

    let mut sample_names = RMQ::new(sample_name_lengths, sample_lcp);

    // Step 1.2: check whether the names are unique.
    if name < sample_tree.len() && rank < k {
        let alpha = name;
        sample_names = sort_tree(
            &mut sample_tree,
            &mut sample_indices,
            &mut sample_names,
            bucket_buffer,
            alpha,
            rank * 3, k
        );
        #[allow(clippy::needless_range_loop)]
        for i in 0..sample_count {
            let index_in_sample_tree = sample_indices[i];
            sample_tree[index_in_sample_tree].name = i + 1;
        }
    }

    // Step 2: sort residual indices.
    //
    // Radix sort using the new names from the sample tree.
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

    macro_rules! write_lcp_and_name_length {
        () => {
            if write_i == 0 {
                output_lcp_values_buffer.push(0);
                output_name_lengths_buffer.push(0);
            } else {
                let previous_node = order[write_i - 1];
                let current_node = order[write_i];
                let (lcp, name_length) = calculate_lcp_and_name_length(
                    tree,
                    &tree_to_sample_tree_mapping,
                    &sample_tree,
                    input_names,
                    &sample_names,
                    residual_modulo,
                    previous_node,
                    current_node
                );
                output_lcp_values_buffer.push(lcp);
                output_name_lengths_buffer.push(name_length);
            }
        };
    }

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

        write_lcp_and_name_length!();
        write_i += 1;
    }
    while residual_i < residual_count {
        let current_residual_node = residual_indices[residual_i];
        order[write_i] = current_residual_node;

        write_lcp_and_name_length!();
        residual_i += 1;
        write_i += 1;
    }
    while sample_i < sample_count {
        let index_in_sample_tree  = sample_indices[sample_i];
        let current_sample_node   = sample_tree[index_in_sample_tree].id;
        order[write_i] = current_sample_node;

        write_lcp_and_name_length!();
        sample_i += 1;
        write_i += 1;
    }

    RMQ::new(output_name_lengths_buffer, output_lcp_values_buffer)
}

#[allow(clippy::too_many_arguments)]
fn calculate_lcp_and_name_length(
    tree: &[Node],
    tree_to_sample_tree_mapping: &[usize],
    sample_tree: &[Node],
    input_names: &RMQ,
    sample_names: &RMQ,
    residual_modulo: usize,
    previous_node: usize,
    current_node: usize,
) -> (usize, usize) {

    #[inline(always)]
    fn lcp_from_tree(tree: &[Node], rmq: &RMQ, previous_node: usize, current_node: usize) -> usize {
        let previous_name = tree[previous_node].name - 1;
        let current_name = tree[current_node].name - 1;
        if previous_name == current_name {
            rmq.name_lengths()[previous_name]
        } else {
            rmq.query(previous_name + 1, current_name)
        }
    }

    let mut lcp = 0;
    let previous_node_modulo = tree[previous_node].depth % 3;
    let current_node_modulo = tree[current_node].depth % 3;

    let previous_name = tree[previous_node].name - 1;
    let current_name = tree[current_node].name - 1;
    let previous_parent = tree[previous_node].parent;
    let current_parent  = tree[current_node].parent;
    let current_name_length = input_names.name_lengths()[current_name];

    let name_length = if current_node_modulo == residual_modulo {
        let current_parent = tree[current_node].parent;
        let current_parent_sample_index = tree_to_sample_tree_mapping[current_parent];
        let current_parent_sample_name = sample_tree[current_parent_sample_index].name - 1;
        current_name_length + sample_names.name_lengths()[current_parent_sample_name]
    } else {
        let current_sample_index = tree_to_sample_tree_mapping[current_node];
        let current_sample_name = sample_tree[current_sample_index].name - 1;
        sample_names.name_lengths()[current_sample_name]
    };

    if previous_node_modulo != residual_modulo && current_node_modulo != residual_modulo {
        let previous_sample_index = tree_to_sample_tree_mapping[previous_node];
        let current_sample_index  = tree_to_sample_tree_mapping[current_node];
        let lcp = lcp_from_tree(
            sample_tree,
            sample_names,
            previous_sample_index,
            current_sample_index,
        );
        return (lcp, name_length);
    }

    // At least one of the nodes is residual now.

    //
    // In all other cases we will compare the node's input name first, therefore, we can make an
    // early exit if those are different.
    //
    if previous_name != current_name {
        let lcp = input_names.query(previous_name + 1, current_name);
        return (lcp, name_length);
    }

    //
    // If both nodes are residual or one of them has depth such that the parent is a sample node,
    // then both nodes' parents are from the sample and we can calculate the lcp value from the
    // sample lcp values.
    //
    if (previous_node_modulo == residual_modulo && current_node_modulo == residual_modulo)
        || ((previous_node_modulo + 1) % 3 == residual_modulo || (current_node_modulo + 1) % 3 == residual_modulo)
    {
        let previous_parent_sample_index = tree_to_sample_tree_mapping[previous_parent];
        let current_parent_sample_index  = tree_to_sample_tree_mapping[current_parent];
        let lcp = current_name_length + lcp_from_tree(
            sample_tree,
            sample_names,
            previous_parent_sample_index,
            current_parent_sample_index,
        );
        return (lcp, name_length);
    }

    //
    // Otherwise, one of the nodes is residual and the other has depth such that the first parent
    // is not a sample node, but the second parent is a sample node.
    //
    let previous_parent_name = tree[previous_parent].name - 1;
    let current_parent_name  = tree[current_parent].name - 1;
    let parent_name_length = input_names.name_lengths()[current_parent_name];
    if previous_parent_name != current_parent_name {
        let lcp = current_name_length + input_names.query(previous_parent_name + 1, current_parent_name);
        return (lcp, name_length);
    }

    let previous_parent_parent              = tree[previous_parent].parent;
    let previous_parent_parent_sample_index = tree_to_sample_tree_mapping[previous_parent_parent];
    let current_parent_parent               = tree[current_parent].parent;
    let current_parent_parent_sample_index  = tree_to_sample_tree_mapping[current_parent_parent];

    let lcp = current_name_length + parent_name_length + lcp_from_tree(
        sample_tree,
        sample_names,
        previous_parent_parent_sample_index,
        current_parent_parent_sample_index,
    );
    (lcp, name_length)
}

#[cfg(test)]
mod tests {
    // TODO: write a test for the decomposition of the SBWT tree to make sure that the order the
    // nodes are traversed is compatible with the level ancestor query support construction
    // procedure.
}

