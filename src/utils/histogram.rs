//! A simple histogram class that can only be incremented.
#![allow(dead_code)]
use anyhow::bail;
use serde::{Deserialize, Serialize};

/// A simple histogram class that can only be incremented. Bins are considered
/// to be 1-based ([0, 1, 2, 3, 4, etc]). Currently, the histogram only supports
/// starting at zero because that is all this package needs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Histogram {
    // Vec-backed value store for the histogram.
    values: Vec<usize>,
    // Starting range for the histogram.
    range_start: usize,
    // Ending range for the histogram.
    range_stop: usize,
}

/// An error that occurs if we try to increment a bin of the histogram that is
/// out-of-bounds for that histogram.
#[derive(Debug, PartialEq, Eq)]
pub struct BinOutOfBoundsError;

impl Histogram {
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
    pub fn range_len(&self) -> usize {
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

    /// Computes the value of the nth percentile based on an exhaustive search.
    pub fn percentile(&self, percentile: f64) -> anyhow::Result<Option<f64>> {
        // (1) Bounds check on the input data
        if !(0.0..=1.0).contains(&percentile) {
            bail!("Provided percentile was not within a valid range.");
        }

        // (2) Count up the total number of items in the histogram
        let mut num_items = 0usize;
        for i in self.range_start..=self.range_stop {
            num_items += self.get(i);
        }

        // (3) If the number of items is zero, then there is no median. Note
        // that num_items is a usize, so obviously cannot be less than zero
        if num_items == 0 {
            return Ok(None);
        }

        // (4) Some simple math to figure out how many items constitutes
        // the nth percentile.
        let needed_items = percentile * num_items as f64;

        // (5) Simple algorithm to calculate the percentile: starting at the
        // lowest value for the histogram, slowly step through the histogram
        // until we have collected `needed_items`.
        let mut collected_items = 0.0;
        let mut index = self.range_start;

        loop {
            // (5a) If we go past the range stop, something really strange has
            // happened and needs to be looked into
            if index > self.range_stop {
                bail!("Unknown error!");
            }

            // (5b) Increment the collected items amount by the current bin we
            // are looking at
            collected_items += self.get(index) as f64;

            // (5c) If the number of collected items eclipses the number of
            // needed items, then we've found our answer!
            if collected_items > needed_items {
                return Ok(Some(index as f64));
            }

            // (5d) If the number of collected items equals the number of needed
            // items, then we have a runoff! Technically, the right way to
            // handle this is to find the next item which has a nonzero count
            // and take the middle of the two (even though that doesn't appear
            // in the set necessarily). So that's what we do here!
            if collected_items == needed_items {
                let lowest = index;
                index += 1;

                while self.get(index) == 0 {
                    index += 1;
                }

                let highest = index;
                return Ok(Some(lowest as f64 + ((highest - lowest) as f64 / 2.0)));
            }

            // (5e) Increment the bin we are looking at by one
            index += 1;
        }
    }

    /// Computes the first quartile of the distribution.
    pub fn first_quartile(&self) -> Option<f64> {
        self.percentile(0.25).unwrap()
    }

    /// Computes the mean of the distribution.
    pub fn median(&self) -> Option<f64> {
        self.percentile(0.5).unwrap()
    }

    /// Computes the third quartile of the distribution.
    pub fn third_quartile(&self) -> Option<f64> {
        self.percentile(0.75).unwrap()
    }

    /// Computes the interquartile range for this distribution.
    pub fn interquartile_range(&self) -> Option<f64> {
        if let Some(first) = self.first_quartile() {
            if let Some(third) = self.third_quartile() {
                return Some(third - first);
            }
        }

        None
    }

    /// Computes the sum of the values within the distribution.
    pub fn sum(&self) -> usize {
        self.values.iter().sum()
    }

    /// Simply returns the values in the distribution by ref.
    pub fn values(&self) -> &[usize] {
        self.values.as_ref()
    }

    /// Normalizes the values so that they sum to one and returns them as a Vec.
    pub fn values_normalized(&self) -> Vec<f64> {
        let total = self.sum() as f64;
        self.values.iter().map(|x| *x as f64 / total).collect()
    }
}

impl Default for Histogram {
    fn default() -> Self {
        Self::zero_based_with_capacity(512)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    pub fn test_initialization() {
        let s = Histogram::zero_based_with_capacity(100);
        assert_eq!(s.range_len(), 101);
        assert_eq!(s.get_range_start(), 0);
        assert_eq!(s.get_range_stop(), 100);
    }

    #[test]
    pub fn test_valid_incremements_and_mean_median() {
        let mut s = Histogram::zero_based_with_capacity(100);
        s.increment(25).unwrap();
        s.increment(50).unwrap();
        s.increment_by(75, 3).unwrap();
        s.increment_by(100, 5).unwrap();

        assert_eq!(s.get(25), 1);
        assert_eq!(s.get(50), 1);
        assert_eq!(s.get(75), 3);
        assert_eq!(s.get(100), 5);

        assert_eq!(s.mean(), 80.0);
        assert_eq!(s.first_quartile().unwrap(), 75.0);
        assert_eq!(s.median().unwrap(), 87.5);
        assert_eq!(s.third_quartile().unwrap(), 100.0);
        assert_eq!(s.interquartile_range().unwrap(), 25.0);
    }

    #[test]
    pub fn test_median_on_empty_histogram() {
        let s = Histogram::zero_based_with_capacity(5000);
        assert!(s.median().is_none());
    }

    #[test]
    pub fn test_median_extensively() {
        let mut s = Histogram::zero_based_with_capacity(5000);

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
        let mut s = Histogram::zero_based_with_capacity(100);
        assert_eq!(s.increment(101).unwrap_err(), BinOutOfBoundsError);
    }

    #[test]
    pub fn test_default_is_one_hundred_zero_based() {
        // The "general" QC facet initializes a SimpleHistogram::default() for
        // its metric collection. As such, it depends on the defaults being
        // from 0 to 100. If you wish to change the default, be sure to update
        // that QC facet accordingly.
        let default = Histogram::default();
        assert_eq!(default.get_range_start(), 0);
        assert_eq!(default.get_range_stop(), 512);
        assert_eq!(default.range_len(), 513);
    }

    #[test]
    pub fn test_values() {
        let mut histogram = Histogram::zero_based_with_capacity(3);
        histogram.increment(1).unwrap();
        histogram.increment(2).unwrap();
        histogram.increment_by(3, 3).unwrap();
        assert_eq!(histogram.values(), [0, 1, 1, 3]);
    }

    #[test]
    pub fn test_values_normalized() {
        let mut histogram = Histogram::zero_based_with_capacity(3);
        histogram.increment(1).unwrap();
        histogram.increment(2).unwrap();
        histogram.increment_by(3, 3).unwrap();
        assert_eq!(histogram.values_normalized(), [0.0, 0.2, 0.2, 0.6]);
    }
}
