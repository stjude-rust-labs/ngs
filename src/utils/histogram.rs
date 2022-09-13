//! A simple histogram class that can only be incremented. Bins are considered
//! to be 1-based ([0, 1, 2, 3, 4, etc]). Currently, the histogram only supports
//! starting at zero because that is all this package needs.
#![allow(dead_code)]
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct SimpleHistogram {
    // Vec-backed value store for the histogram.
    values: Vec<usize>,
    // Starting range for the histogram.
    range_start: usize,
    // Ending range for the histogram.
    range_stop: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub struct BinOutOfBoundsError;

impl SimpleHistogram {
    /// Creates a zero-based histogram with a given capacity.
    pub fn zero_based_with_capacity(capacity: usize) -> Self {
        Self {
            values: vec![0; capacity + 1],
            range_start: 0,
            range_stop: capacity,
        }
    }

    /// Increments a particular bin in the histogram by the specified value.
    pub fn increment_by(&mut self, bin: usize, value: usize) -> Result<(), BinOutOfBoundsError> {
        if bin < self.range_start || bin > self.range_stop {
            return Err(BinOutOfBoundsError);
        }

        self.values[bin] += value;
        Ok(())
    }

    /// Increments a particular bin in the histogram by one.
    pub fn increment(&mut self, bin: usize) -> Result<(), BinOutOfBoundsError> {
        self.increment_by(bin, 1)
    }

    /// Gives the length of the open range for the histogram.
    pub fn len(&self) -> usize {
        self.range_stop - self.range_start + 1
    }

    /// Gives the starting position for the open range of the histogram.
    pub fn get_range_start(&self) -> usize {
        self.range_start
    }

    /// Gives the stopping position for the open range of the histogram.
    pub fn get_range_stop(&self) -> usize {
        self.range_stop
    }

    /// Gets a value for a bin within a histogram.
    pub fn get(&self, bin: usize) -> usize {
        *self
            .values
            .get(bin)
            .expect("Could not lookup value for template length histogram bin.")
    }

    /// Computes the mean of all values within the histogram.
    pub fn mean(&self) -> f64 {
        let mut sum = 0.0;
        let mut denominator = 0.0;

        for i in self.range_start..=self.range_stop {
            let bin_value = self.get(i);
            denominator += bin_value as f64;
            sum += (bin_value * i) as f64;
        }

        sum / denominator
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    pub fn test_initialization() {
        let s = SimpleHistogram::zero_based_with_capacity(100);
        assert_eq!(s.len(), 101);
        assert_eq!(s.get_range_start(), 0);
        assert_eq!(s.get_range_stop(), 100);
    }

    #[test]
    pub fn test_valid_incremements_and_mean() {
        let mut s = SimpleHistogram::zero_based_with_capacity(100);
        s.increment(25).unwrap();
        s.increment(50).unwrap();
        s.increment_by(75, 3).unwrap();

        assert_eq!(s.get(25), 1);
        assert_eq!(s.get(50), 1);
        assert_eq!(s.get(75), 3);

        assert_eq!(s.mean(), 60.0);
    }

    #[test]
    pub fn test_invalid_increments() {
        let mut s = SimpleHistogram::zero_based_with_capacity(100);
        assert_eq!(s.increment(101).unwrap_err(), BinOutOfBoundsError);
    }
}
