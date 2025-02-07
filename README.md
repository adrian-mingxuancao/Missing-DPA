# Missing DPA Analysis

This repository contains implementations and comparisons of different DPA (Deterministic Population Analysis) methods for handling missing data.

## Methods

The project implements and compares two main methods:

1. **DPA-MAR-3**: Original method using sequential testing with empirically adjusted thresholds
2. **DPA-MAR-Limit**: New method using limiting distribution theory

### Key Differences

The main differences between DPA-MAR-Limit and DPA-MAR-3 are:

1. **Threshold Calculation**: 
   - DPA-MAR-Limit uses `(1/p0)*(1 + sqrt(gamma))^2 - ((1-p0)/p0)`
   - DPA-MAR-3 uses the `upper_edge` function with empirical adjustments

2. **Missing Value Handling**:
   - DPA-MAR-Limit explicitly calculates observed probability p0 from data
   - DPA-MAR-3 uses direct adjustment based on input missing probability

3. **Scaling**:
   - DPA-MAR-Limit implements custom scaling for missing values
   - DPA-MAR-3 uses R's standard `scale()` function

## Project Structure

- `Code/`: Contains R implementation files
  - `Mar_Rank.R`: Core functions for rank estimation
  - `dpa_mar_comparison.R`: Comparison script and implementations
  - `analyze_geo_data.R`: Analysis of GEO dataset
- `experiment_data/`: Contains experimental results

## Usage

To run the comparison between methods:

```R
source("Code/Mar_Rank.R")
source("Code/dpa_mar_comparison.R")
```

## Results

Our experiments show that:
1. DPA-MAR-3 is more robust to high missing probabilities
2. DPA-MAR-Limit performs better with low missing probabilities and high signal strength
3. Both methods perform well for lower ranks (3-5) but diverge in performance for higher ranks (8)
