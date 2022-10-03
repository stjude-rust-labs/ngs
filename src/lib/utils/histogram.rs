//! A simple histogram class that can only be incremented. Bins are considered
//! to be 1-based ([0, 1, 2, 3, 4, etc]). Currently, the histogram only supports
//! starting at zero because that is all this package needs.
#![allow(dead_code)]
use anyhow::bail;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
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

    /// Gives the stopping position for the open range of the histogram.
    pub fn in_range(&self, value: usize) -> bool {
        self.range_start <= value && value <= self.range_stop
    }

    /// Gets a value for a bin within a histogram.
    pub fn get(&self, bin: usize) -> usize {
        *self.values.get(bin).unwrap_or_else(|| {
            panic!(
                "Could not lookup value for template length histogram bin: {}.",
                bin
            )
        })
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

    /// Computes the median of all values within the histogram.
    pub fn median(&self) -> Option<f64> {
        let mut sum: i64 = 0;
        // fp => Front pointer
        // bp => Back pointer
        let mut fp = self.range_start as i64 - 1;
        let mut bp = self.range_stop as i64 + 1;
        let mut last_known_nonzero_front: Option<i64> = None;
        let mut last_known_nonzero_back: Option<i64> = None;

        while fp != bp {
            if sum < 0 {
                fp += 1;

                if !self.in_range(fp as usize) {
                    break;
                }

                let val = self.get(fp as usize);
                if val != 0 && fp != bp {
                    last_known_nonzero_front = Some(fp);
                    sum += val as i64;
                }
            } else {
                bp -= 1;

                if !self.in_range(bp as usize) {
                    break;
                }

                let val = self.get(bp as usize);
                if val != 0 && fp != bp {
                    last_known_nonzero_back = Some(bp);
                    sum -= val as i64;
                }
            }
        }

        if sum == 0 {
            if let Some(nonzero_front) = last_known_nonzero_front {
                if let Some(nonzero_back) = last_known_nonzero_back {
                    return Some(
                        (nonzero_back - nonzero_front) as f64 / 2.0 + nonzero_front as f64,
                    );
                }
            }

            None
        } else {
            Some(fp as f64)
        }
    }
}

impl Default for SimpleHistogram {
    fn default() -> Self {
        Self::zero_based_with_capacity(512)
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
    pub fn test_valid_incremements_and_mean_median() {
        let mut s = SimpleHistogram::zero_based_with_capacity(100);
        s.increment(25).unwrap();
        s.increment(50).unwrap();
        s.increment_by(75, 3).unwrap();

        assert_eq!(s.get(25), 1);
        assert_eq!(s.get(50), 1);
        assert_eq!(s.get(75), 3);

        assert_eq!(s.mean(), 60.0);
        assert_eq!(s.median().unwrap(), 75.0);
    }

    #[test]
    pub fn test_median_on_empty_histogram() {
        let s = SimpleHistogram::zero_based_with_capacity(5000);
        assert!(s.median().is_none());
    }

    #[test]
    pub fn test_median_extensively() {
        let mut s = SimpleHistogram::zero_based_with_capacity(5000);

        // Start to add in values
        s.increment_by(0, 2500).unwrap();
        s.increment_by(10, 2500).unwrap();
        s.increment_by(100, 2500).unwrap();
        s.increment_by(5000, 5000).unwrap();
        let median = s.median();
        assert!(median.is_some());
        assert_eq!(median.unwrap(), 100.0);

        // If there is a tie, take the value in between the two middle values
        s.increment_by(200, 2500).unwrap();
        let median = s.median();
        assert!(median.is_some());
        assert_eq!(median.unwrap(), 150.0);

        // If we add one more to sway the vote, should shift the median
        s.increment(200).unwrap();
        let median = s.median();
        assert!(median.is_some());
        assert_eq!(median.unwrap(), 200.0);
    }

    #[test]
    pub fn test_invalid_increments() {
        let mut s = SimpleHistogram::zero_based_with_capacity(100);
        assert_eq!(s.increment(101).unwrap_err(), BinOutOfBoundsError);
    }

    #[test]
    pub fn test_default_is_one_hundred_zero_based() {
        // The "general" QC facet initializes a SimpleHistogram::default() for
        // its metric collection. As such, it depends on the defaults being
        // from 0 to 100. If you wish to change the default, be sure to update
        // that QC facet accordingly.
        let default = SimpleHistogram::default();
        assert_eq!(default.get_range_start(), 0);
        assert_eq!(default.get_range_stop(), 512);
        assert_eq!(default.len(), 513);
    }
}
