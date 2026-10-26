// https://codeforces.com/blog/entry/78931

pub const BLOCK_SIZE: usize = usize::BITS as usize;
pub struct RangeMinimumSupport {
    values: Vec<usize>,
    masks: Vec<usize>,
    sparse_table: Vec<usize>,
}

impl RangeMinimumSupport {
    pub fn new(values: Vec<usize>) -> RangeMinimumSupport {
        let length = values.len();
        let mut masks  = Vec::<usize>::with_capacity(length);
        let mut sparse_table = vec![0_usize; length];

        let mut current_mask = 0_usize;

        for i in 0..length {
            current_mask = current_mask.unbounded_shl(1);

            while current_mask > 0 {
                let offset = current_mask.trailing_zeros() as usize;
                if i == Self::_index_of_smaller(&values, i, i - offset) {
                    current_mask ^= 1 << offset;
                } else {
                    break;
                }
            }

            current_mask |= 1;
            masks.push(current_mask);
        }

        let block_count = length / BLOCK_SIZE;
        #[allow(clippy::needless_range_loop)]
        for column in 0..block_count {
            let last_position_in_block = column * BLOCK_SIZE + BLOCK_SIZE - 1;
            sparse_table[column] = Self::_in_block_index_query(
                &masks,
                last_position_in_block,
                BLOCK_SIZE
            );
        }
        let mut row = 1;
        while (1_usize << row) <= block_count {
            let current_row_offset = block_count * row;
            let previous_row_offset = current_row_offset - block_count;
            let last_column = block_count - (1 << row);
            let neighbour_block_offset = 1 << (row - 1);
            for column in 0..=last_column {
                sparse_table[current_row_offset + column] = Self::_index_of_smaller(
                    &values,
                    sparse_table[previous_row_offset + column],
                    sparse_table[previous_row_offset + column + neighbour_block_offset]
                );
            }
            row += 1;
        }

        RangeMinimumSupport {
            values,
            sparse_table,
            masks,
        }
    }

    fn _index_of_smaller(values: &[usize], a: usize, b: usize) -> usize {
        if values[a] <= values[b] { a } else { b }
    }

    fn index_of_smaller(&self, a: usize, b: usize) -> usize {
        Self::_index_of_smaller(&self.values, a, b)
    }

    fn _in_block_index_query(masks: &[usize], right: usize, length: usize) -> usize {
        let mask: usize = 1_usize.unbounded_shl(length as u32).wrapping_sub(1);
        let offset = BLOCK_SIZE - 1 - (masks[right] & mask).leading_zeros() as usize;
        right - offset
    }

    fn in_block_index_query(&self, right: usize, length: usize) -> usize {
        Self::_in_block_index_query(&self.masks, right, length)
    }

    pub fn query(&self, l: usize, r: usize) -> usize {
        #[allow(clippy::int_plus_one)]
        if r - l + 1 <= BLOCK_SIZE {
            let index = self.in_block_index_query(r, r - l + 1);
            return self.values[index];
        }

        let mut answer = Self::_index_of_smaller(
            &self.values,
            self.in_block_index_query(l + BLOCK_SIZE - 1, BLOCK_SIZE),
            self.in_block_index_query(                 r, BLOCK_SIZE),
        );

        let left_block  = l / BLOCK_SIZE + 1;
        let right_block = r / BLOCK_SIZE - 1;

        if left_block <= right_block {
            let block_range_length = right_block - left_block + 1;
            let row = BLOCK_SIZE - 1 - block_range_length.leading_zeros() as usize;
            let row_block_range_length = (1 << row) - 1;
            let block_count = self.values.len() / BLOCK_SIZE;
            answer = self.index_of_smaller(
                answer,
                self.index_of_smaller(
                    self.sparse_table[row * block_count + left_block],
                    self.sparse_table[row * block_count + right_block - row_block_range_length],
                )
            );
        }

        self.values[answer]
    }

    pub fn values(&self) -> &[usize] {
        &self.values
    }

    pub fn swap_values(&mut self, mut values: Vec<usize>) -> Vec<usize> {
        std::mem::swap(&mut self.values, &mut values);
        values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rmq() {
        let values = (1..=512).collect::<Vec<usize>>();
        let rms = RangeMinimumSupport::new(values.clone());
        for i in 0..512 {
            for j in i..512 {
                let mut real_min_value = 1000;
                #[allow(clippy::needless_range_loop)]
                for k in i..=j {
                    real_min_value = real_min_value.min(values[k]);
                }
                let query_min_value = rms.query(i, j);
                assert_eq!(real_min_value, query_min_value, "{}-{}", i, j);
            }
        }
    }
}

