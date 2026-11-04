// https://codeforces.com/blog/entry/126580
// https://judge.yosupo.jp/submission/242808

use super::Node;

#[derive(Debug)]
struct NodeInfo {
    ladder_start: usize,
    ladder_length: usize,
    euler_index: usize,
    deepest_leaf: usize,
}

#[derive(Debug)]
struct LevelAncestorSupport {
    jump_pointers: Vec<(usize, usize)>,
    ladder_values: Vec<usize>,
    auxiliary_info: Vec<NodeInfo>,
}

impl LevelAncestorSupport {
    fn new(tree: &[Node]) -> Self {
        let tree_size = tree.len();

        let mut jump_pointers  = Vec::<(usize, usize)>::with_capacity(2 * tree_size);
        let mut ladder_values  = Vec::<usize>::with_capacity(3 * tree_size);
        let mut auxiliary_info = Vec::<NodeInfo>::with_capacity(tree_size);

        let mut stack = Vec::<usize>::new();

        // simulate dfs on the tree and find the deepest leaves
        // populate the jump pointers
        for (node_index, node) in tree.iter().enumerate() {
            while let Some(last) = stack.last() && *last > node.parent {
                // update deepest leaf
                stack.pop();
            }

            stack.push(node_index);

        }
 
        // calculate the ladders

        todo!()
    }

    fn query(node: usize, depth: usize) -> usize {
        todo!()
    }
}

