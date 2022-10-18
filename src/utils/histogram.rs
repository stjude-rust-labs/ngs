//! Histogram used as the basis for statistics counting in many quality control
//! facets.
//!
//! # Overview
//!
//! [Histograms] are a common way to make sense of a distribution of data.
//! Briefly, data along a continuous (or sometimes discrete) distribution is
//! partitioned into bins of a predetermined size. Bins are generally
//! consecutive and are non-overlapping in nature. This straightforward model
//! helps to easily visualize and make sense of sometimes complex data.
//!
//! In this module, we implement a very simple histogram that represents the
//! minimum viable histogram needed for the `ngs` command line tool. Commonly,
//! this is to count discrete data that falls within a particular numerical
//! distribution. Simply put, the histogram follows these rules:
//!
//! 1. Only discrete numbers are considered as bins. In other words, bins
//!    represent values in the range of `[0, 1, 2, 3, ..., n]`.
//! 2. The range of numerical values always starts at zero and ends at the
//!    capacity specified in the constructor. This is because all of the metrics
//!    we have collected within a histogram (to this point) should include zero
//!    as a valid bin.
//!
//! As we continue to develop the `ngs` command line tool, these characteristics
//! may change—please treat these contracts as transient!
//!
//! # Usage
//!
//! You can create a histogram with only one initialization function:
//! [`zero_based_with_capacity`][Histogram::zero_based_with_capacity]. This
//! design is intentional for now, we only want you creating histograms that
//! start at zero!
//!
//! After you have a [`Histogram`], you can start doing things with it!
//! Commonly, you'll want to increment bins within a histogram: you can
//! increment a bin by one (essentially, a `+= 1`) with the
//! [`increment`][`Histogram::increment`] method. If you need to increment the
//! bin by more than one at a time, you can use the
//! [`increment_by`][Histogram::increment_by] method.
//!
//! ```
//! use ngs::utils::histogram::Histogram;
//! let mut hist = Histogram::zero_based_with_capacity(10);
//!
//! // Increments the zero bin by one.
//! let result = hist.increment(0);
//! assert!(result.is_ok());
//!
//! // Increments the one bin by fourty-two.
//! let result = hist.increment_by(1, 42);
//! assert!(result.is_ok());
//!
//! // Ensure that we actually recorded these values.
//! assert_eq!(hist.get(0), 1);
//! assert_eq!(hist.get(1), 42);
//! ```
//!
//! If necessary, you can examine the range of the histogram:
//!
//! ```
//! use ngs::utils::histogram::Histogram;
//! let mut hist = Histogram::zero_based_with_capacity(10);
//!
//! // Remember, this includes zero!
//! assert_eq!(hist.range_len(), 11);
//! assert_eq!(hist.range_start(), 0);
//! assert_eq!(hist.range_stop(), 10);
//! ```
//!
//! You can even check if a value falls within the range of the [`Histogram`]!
//!
//! ```
//! use ngs::utils::histogram::Histogram;
//! let mut hist = Histogram::zero_based_with_capacity(10);
//!
//! assert!(hist.in_range(0));
//! assert!(hist.in_range(5));
//! assert!(hist.in_range(10));
//! assert!(!hist.in_range(11));
//! ```
//!
//! Note that if you try to increment the [`Histogram`] within a bin that falls
//! outside its range, you will get a [`BinOutOfBoundsError`].
//!
//! ```
//! use ngs::utils::histogram::Histogram;
//! use ngs::utils::histogram::BinOutOfBoundsError;
//! let mut hist = Histogram::zero_based_with_capacity(10);
//!
//! let result = hist.increment(11);
//! assert_eq!(result.unwrap_err(), BinOutOfBoundsError);
//! ```
//!
//! Finally, if all you want are the values out of the histogram (whether raw or
//! normalized), you can use the [`values`][Histogram::values] and
//! [`values_normalized`][Histogram::values_normalized] methods respectively:
//!
//! ```
//! use ngs::utils::histogram::Histogram;
//! let mut hist = Histogram::zero_based_with_capacity(3);
//!
//! hist.increment_by(0, 2).unwrap();
//! hist.increment(1).unwrap();
//! hist.increment(2).unwrap();
//!
//! assert_eq!(hist.values(), [2, 1, 1, 0]);
//! assert_eq!(hist.values_normalized(), [0.5, 0.25, 0.25, 0.0]);
//! ```
//!
//! You can do other various operations, such as:
//!
//! - Find the mean of the distribution ([`mean`][Histogram::mean]).
//! - Find an arbitrary percentile of the distribution ([`percentile`][Histogram::percentile]).
//! - Find the first quartile of the distribution ([`first_quartile`][Histogram::first_quartile]).
//! - Find the median (second quartile) of the distribution ([`median`][Histogram::median]).
//! - Find the third quartile of the distribution ([`third_quartile`][Histogram::third_quartile]).
//! - Find the interquartile range of the distribution ([`interquartile_range`][Histogram::interquartile_range]).
//! - Find the sum of all counts within the distribution ([`sum`][Histogram::sum]).
//!
//! [Histograms]: https://en.wikipedia.org/wiki/Histogram
//!
//! ## Why not the pre-existing `histogram` crate?
//!
//! Some consideration was given to using the (seemingly) quite good
//! [`histogram`] crate that is already published on [crates.io]. When it boiled
//! down to it, the authors choices were:
//!
//! 1. Go and make a PR to add `Serialize` and `Deserialize` to the existing
//!    crate.
//! 2. Create our own customized histogram, or
//!
//! The first option had no guarantee of working out—especially since the last
//! commit date was over two years ago from the time of writing. Further, the
//! layout of the struct they used did not seem to be optimal for our
//! serialization/deserialization needs. Thus, although that crate has many
//! methods that are similar (and sometimes inspired our approach within this
//! module), we ultimately decided just to roll our own.
//!
//! [`histogram`]: https://docs.rs/histogram/latest/histogram/
//! [crates.io]: https://crates.io

use anyhow::bail;
use serde::{Deserialize, Serialize};

/// Histogram used as the basis for statistics counting in many quality control
/// facets. For more in depth information, please see the [module-level
/// documentation].
///
/// [module-level documentation]: self
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
    //=================//
    // Initializations //
    //=================//

    /// Creates a zero-based histogram with a given capacity.
    pub fn zero_based_with_capacity(capacity: usize) -> Self {
        Self {
            values: vec![0; capacity + 1],
            range_start: 0,
            range_stop: capacity,
        }
    }

    //=================================//
    // Getting and incrementing values //
    //=================================//

    /// Increments a particular bin in the histogram by one.
    pub fn increment(&mut self, bin: usize) -> Result<(), BinOutOfBoundsError> {
        self.increment_by(bin, 1)
    }

    /// Increments a particular bin in the histogram by the specified value.
    pub fn increment_by(&mut self, bin: usize, value: usize) -> Result<(), BinOutOfBoundsError> {
        if bin < self.range_start || bin > self.range_stop {
            return Err(BinOutOfBoundsError);
        }

        self.values[bin] += value;
        Ok(())
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

    /// Simply returns the values in the distribution by ref.
    pub fn values(&self) -> &[usize] {
        self.values.as_ref()
    }

    /// Normalizes the values so that they sum to one and returns them as a Vec.
    pub fn values_normalized(&self) -> Vec<f64> {
        let total = self.sum() as f64;
        self.values.iter().map(|x| *x as f64 / total).collect()
    }

    //=======//
    // Range //
    //=======//

    /// Gives the length of the open range for the histogram.
    pub fn range_len(&self) -> usize {
        self.range_stop - self.range_start + 1
    }

    /// Gives the starting position for the open range of the histogram.
    pub fn range_start(&self) -> usize {
        self.range_start
    }

    /// Gives the stopping position for the open range of the histogram.
    pub fn range_stop(&self) -> usize {
        self.range_stop
    }

    /// Indicates whether a particular value falls within the range of the histogram.
    ///
    /// ```
    /// use ngs::utils::histogram::Histogram;
    /// let mut hist = Histogram::zero_based_with_capacity(100);
    ///
    /// assert!(hist.in_range(0));
    /// assert!(hist.in_range(100));
    /// assert!(!hist.in_range(101));
    /// ```
    pub fn in_range(&self, value: usize) -> bool {
        (self.range_start..=self.range_stop).contains(&value)
    }

    //========================//
    // Numerical computations //
    //========================//

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

    /// Counts from the bottom of the histogram up until a certain bin.
    pub fn count_from_bottom_until(&self, bin: usize) -> usize {
        let mut sum = 0;
        for i in self.range_start()..=bin {
            sum += self.get(i)
        }

        sum
    }

    /// Counts from the top of the histogram down until a certain bin.
    ///
    /// More accurately, we technically walk from the specified bin up to the
    /// top of the histogram, but this provides the same results.
    pub fn count_from_top_until(&self, bin: usize) -> usize {
        let mut sum = 0;
        for i in bin..=self.range_stop() {
            sum += self.get(i)
        }

        sum
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
        assert_eq!(s.range_start(), 0);
        assert_eq!(s.range_stop(), 100);
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
        assert_eq!(default.range_start(), 0);
        assert_eq!(default.range_stop(), 512);
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

    #[test]
    pub fn test_count_values_from_bottom() {
        let mut histogram = Histogram::zero_based_with_capacity(3);
        histogram.increment_by(0, 5).unwrap();
        histogram.increment_by(1, 3).unwrap();
        histogram.increment_by(2, 6).unwrap();
        assert_eq!(histogram.count_from_bottom_until(0), 5);
        assert_eq!(histogram.count_from_bottom_until(1), 8);
        assert_eq!(histogram.count_from_bottom_until(2), 14);
        assert_eq!(histogram.count_from_bottom_until(3), 14);
    }

    #[test]
    pub fn test_count_values_from_top() {
        let mut histogram = Histogram::zero_based_with_capacity(3);
        histogram.increment_by(0, 5).unwrap();
        histogram.increment_by(1, 3).unwrap();
        histogram.increment_by(2, 6).unwrap();
        assert_eq!(histogram.count_from_top_until(3), 0);
        assert_eq!(histogram.count_from_top_until(2), 6);
        assert_eq!(histogram.count_from_top_until(1), 9);
        assert_eq!(histogram.count_from_top_until(0), 14);
    }
}
