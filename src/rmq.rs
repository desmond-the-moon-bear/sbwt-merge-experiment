
pub struct RangeMinimumSupport {
    pub length: usize,
    capacity: usize,
    blocks: Vec<usize>,
    masks: Vec<usize>,
}

impl RangeMinimumSupport {
    pub fn with_capacity(capacity: usize) -> RangeMinimumSupport {
        RangeMinimumSupport {
            length: 0,
            capacity,
            blocks: vec![0_usize; capacity],
            masks: vec![0_usize; capacity],
        }
    }

    pub fn build(&mut self, input: &[usize]) {
        self.length = input.len();
        // TODO: ...
    }

    pub fn query(&self, l: usize, r: usize) -> usize {
        // TODO: ...
        0
    }
}

