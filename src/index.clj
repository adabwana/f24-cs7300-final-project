^:kindly/hide-code
(ns index
  (:require
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind])
  (:import [java.time LocalDate]
           [java.time.format DateTimeFormatter]))

^:kindly/hide-code
(def md (comp kindly/hide-code kind/md))

(let [formatter (DateTimeFormatter/ofPattern "M/d/yy")
      current-date (str (.format (LocalDate/now) formatter))]
  (md (str "
            
### Jaryt Salvo
**Date:** **" current-date "**

**Fall 2024 | CS 7300 Unsupervised Learning**

*************

# Statistical Computing: A Functional Approach

## Foundations and Implementation

This work presents a rigorous implementation of fundamental *statistical algorithms* through the lens of *functional programming*. By leveraging `Clojure`'s *immutable data structures* and *pure functions*, we develop a robust framework for *numerical computing* that emphasizes both mathematical precision and computational efficiency.

## Descriptive Statistics

Statistical measures form the cornerstone of *numerical computing*. Our implementation establishes fundamental operations with careful attention to *numerical stability* and computational efficiency. The `descriptive` namespace provides:

- Robust implementations of *central tendency* and *dispersion measures*
- Numerically stable algorithms for *variance computation*
- Efficient *matrix operations* for *covariance analysis*

## Eigendecomposition

The `eigen` namespace demonstrates the synthesis of *mathematical theory* with practical computation through:

- **Power iteration** for dominant *eigenvalue computation*
- **QR algorithm** for complete *eigendecomposition*
- **Inverse iteration** for *eigenvector determination*

## Principal Component Analysis

The `pca` namespace builds directly on these foundations, implementing **PCA** with:

- Theoretical development from *covariance analysis*
- Efficient *matrix transformations* for high-dimensional data
- Comprehensive comparison with industry-standard implementations

## Implementation Validation

Systematic validation ensures correctness through comparison with established tools:

- Baseline implementation in `Python` using `scikit-learn`
- Comprehensive test suite ensuring *numerical accuracy*
- Performance analysis and optimization strategies

## Technical Architecture

The implementation leverages modern `Clojure` libraries while maintaining *functional purity*:

- `Neanderthal` for efficient *numerical operations*
- Property-based testing for robust validation
- Comprehensive documentation integrated with code

This work serves dual purposes: as a practical implementation of *statistical algorithms* and as an exploration of *functional programming*'s capabilities in *numerical computing*. Through careful attention to both mathematical rigor and computational efficiency, we demonstrate the viability of *functional programming* for serious statistical computing.")))