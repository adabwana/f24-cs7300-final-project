(ns neandersolve.descriptive
  (:require
   [criterium.core :refer [quick-bench]]
   [uncomplicate.fluokitten.core :refer [fmap!]]
   [uncomplicate.commons.core :refer [let-release with-release]]
   [uncomplicate.neanderthal
    [core :refer [col copy dim dot entry entry! mrows mv!
                  ncols nrm2 scal! sum trans vctr]]
    [native :refer [dge dv fge fv]]
    [math :refer [exp sqr sqrt]]
    [vect-math :refer [abs! linear-frac! log!]]]))

; # Descriptive Statistics

; Statistical analysis begins with understanding the fundamental characteristics of our data. These *descriptive measures* help us summarize and interpret large datasets through key numerical values. Our implementation prioritizes both *numerical stability* and computational efficiency.

; ## Basic Statistical Measures

; ### Arithmetic Mean

; The *arithmetic mean* forms the foundation of many statistical computations. For a vector $x = (x_1, ..., x_n)$, the mean is defined as:

; $$\bar{x} = \frac{1}{n}\sum_{i=1}^n x_i$$

; While simple in theory, implementation requires careful consideration of *numerical stability* and efficiency.

(defn mean
  "Computes the arithmetic mean of a vector using pure functional operations.
   Returns a scalar value representing the mean."
  ^double [x]
  (/ (sum x) (dim x)))

; For weighted data, we extend this concept to the *weighted mean*:

; $$\bar{x}_w = \frac{\sum_{i=1}^n w_i x_i}{\sum_{i=1}^n w_i}$$

; where $w_i$ represents the weight of each observation. This implementation handles edge cases by returning `NaN` for invalid inputs.

(defn weighted-mean
  "Computes the weighted mean of a vector given a weight vector.
   Both vectors must have the same dimension.
   Returns NaN if the sum of weights is zero or if vectors have different dimensions."
  ^double [x weights]
  (let [n (dim x)
        weight-sum (sum weights)]
    (if (and (= n (dim weights))
             (not (zero? weight-sum)))
      (/ (dot x weights) weight-sum)
      Double/NaN)))

; ### Geometric Mean

; The *geometric mean* provides a natural measure of *central tendency* for multiplicative processes and relative changes. For a vector $x$ of positive numbers, it is defined as:

; $$\bar{x}_g = \left(\prod_{i=1}^n x_i\right)^{1/n} = \exp\left(\frac{1}{n}\sum_{i=1}^n \ln x_i\right)$$

; We compute this by transforming to *log space* to avoid numerical overflow:

(defn geometric-mean!
  "Computes the geometric mean of a vector using pure functional operations.
   Returns a scalar value representing the geometric mean."
  ^double [x!] ; (abs! x!) to allow negatives
  (exp (mean (log! x!))))

; ### Harmonic Mean

; The *harmonic mean* is particularly useful for averaging rates and speeds, giving more weight to smaller values. For a vector $x$, it is defined as:

; $$\bar{x}_h = \frac{n}{\sum_{i=1}^n \frac{1}{x_i}} = n \left( \sum_{i=1}^n \frac{1}{x_i} \right)^{-1}$$

(defn harmonic-mean!
  "Computes the harmonic mean of a vector.
   Returns NaN for vectors containing zeros.
   Modifies the input vector in place."
  ^double [x!]
  (let [n (dim x!)]
    (if (< 0 n)
      (/ n (sum (fmap! #(/ 1.0 %) x!)))
      Double/NaN)))

; The arithmetic (AM), geometric (GM), and harmonic means (HM) are related through a fundamental inequality:

; $$ \mathrm{AM} \geq \mathrm{GM} \geq \mathrm{HM} $$

; This relationship, known as the AM-GM-HM inequality, demonstrates a consistent ordering among these measures of central tendency. The equality case occurs if and only if all elements in the sample are identical, reflecting the means' convergence for homogeneous data. See [Wikipedia](https://en.wikipedia.org/wiki/Mean#Relationship_between_AM,_GM,_and_HM).

; ## Variance Computation

; The *sample variance* measures the spread of data around its mean. While the theoretical formula is straightforward:

; $$\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2$$

; A numerically stable implementation requires careful consideration to avoid *catastrophic cancellation*. We use a two-pass algorithm that first centers the data:

(defn variance!
  "Computes the sample variance using a numerically stable algorithm.
   Modifies the input vector x! in place.
   Returns NaN for empty vectors."
  [x!]
  (let [n (dim x!)]
    (if (< 0 n)
      (/ (sqr (nrm2 (linear-frac! 1.0 x! (- (mean x!))))) (dec n))
      Double/NaN)))

; Our implementation computes *variance* using vector operations for efficiency:

; - `linear-frac!` centers the data by subtracting the mean: $(x_i - x̄)$
; - `nrm2` computes the *Euclidean norm*: $\sqrt{\sum (x_i - x̄)^2}$
; - `sqr` gives us the sum of squared deviations: $\sum (x_i - x̄)^2$
; - Finally, we divide by $(n-1)$ for the unbiased *sample variance*

; This *vectorized approach* is mathematically equivalent to:
; $$\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2$$


; For cases where *preserving* the input vector is important, we provide a *non-destructive version* that works on a copy of the data:

(defn variance
  "Computes the sample variance using a numerically stable algorithm.
   Preserves the input vector by working on a copy."
  [x]
  (with-release [x-copy (copy x)]
    (variance! x-copy)))

(variance (fv 1 2 3 4 5 6 7 8 9 10))

; ## Mean and Variance in a Single Pass

; While separate computations of mean and variance are straightforward, we can improve efficiency by computing both statistics in a single pass through the data. This approach reduces memory access and computational overhead, particularly important for large datasets.

; The algorithm maintains running sums for both the first and second *moments* of the data. For a dataset of size $n$, we compute:

; $$\mu = \frac{1}{n}\sum_{i=1}^n x_i$$
; $$\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (x_i - \mu)^2$$

; The implementation modifies the input vector in place for efficiency, returning both statistics as a pair:

(defn mean-variance!
  "Computes both mean and variance in a single pass.
   Modifies the input vector x! in place.
   Returns [mean variance] pair, or [NaN NaN] for empty vectors."
  [x!]
  (let [n (dim x!)]
    (if (< 0 n)
      (let [mu (mean x!)]
        [mu (/ (sqr (nrm2 (linear-frac! 1.0 x! (- mu)))) (dec n))]) ; or n?
      [Double/NaN Double/NaN])))

(mean-variance! (fv 1 2 3 4 5 6 7 8 9 10))

; ## Matrix Operations

; Statistical computations often require working with matrices, where each column represents a variable and each row an observation. We begin with computing means across *matrix columns*.

; The *matrix mean* operation computes the *arithmetic mean* for each column, effectively reducing a matrix to a vector of column means:

(defn matrix-mean
  "Computes column means of a matrix.
   Returns a vector where each element is the mean of the corresponding column."
  ([a ones res!]
   (if (< 0 (mrows a))
     (mv! (/ 1.0 (mrows a)) (trans a) ones 0.0 res!)
     (entry! res! Double/NaN)))
  ([a ones]
   (let-release [res (vctr a (ncols a))]
                (matrix-mean a ones res)))
  ([a]
   (with-release [ones (entry! (vctr a (mrows a)) 1.0)]
     (matrix-mean a ones))))

(matrix-mean (fge 3 2 [4 2 3 4 5 6]))

; The matrix variance computation extends our vector operations to handle multiple variables simultaneously. For each column $j$ in matrix $A$, we compute:

; $$\sigma^2_j = \frac{1}{n-1}\sum_{i=1}^n (a_{ij} - \mu_j)^2$$

; We provide both in-place and copying versions to accommodate different usage patterns:

(defn matrix-variance!
  "Computes column variances of a matrix in-place.
   Modifies the input matrix a! during computation."
  [a!]
  (let [n (ncols a!)
        result (vctr a! n)]
    (loop [j 0]
      (when (< j n)
        (entry! result j (variance (col a! j)))
        (recur (inc j))))
    result))

(defn matrix-variance
  "Computes column variances of a matrix.
   Preserves the input matrix by working on a copy."
  [a]
  (with-release [a-copy (copy a)]
    (matrix-variance! a-copy)))

; ## Standard Error

; The standard error of the mean quantifies the uncertainty in our estimate of the population mean. For a sample of size $n$, it is defined as:

; $$SE = \frac{\sigma}{\sqrt{n}}$$

(defn standard-error
  "Computes the standard error of the mean."
  [x]
  (let [n (dim x)]
    (sqrt (/ (variance x) n))))

; It is important to distinguish between the standard error of the mean and the standard deviation. While both are frequently used to summarize data, they serve distinct statistical purposes:

; - The standard deviation characterizes the variability of individual observations within the sample
; - The standard error quantifies the precision of the sample mean as an estimator of the population mean

; As sample size increases, the standard error decreases (approaching zero), reflecting improved estimation precision, while the standard deviation converges to the true population parameter. This relationship emerges from the central limit theorem and forms the basis for statistical inference.

; The standard deviation, $\sigma$, is the sample standard deviation calculated as the square root of the variance:

; $$\sigma = \sqrt{\mathrm{Var}(X)}$$

(defn sd [x]
  (sqrt (variance x)))

; ## Covariance Calculations

; *Covariance* measures the joint variability between two variables. For vectors $x$ and $y$, the *sample covariance* is defined as:

; $$cov(x,y) = \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})$$

; The implementation first centers both vectors and then computes their *dot product*:

(defn covariance!
  "Computes covariance between two vectors in-place.
   Modifies both input vectors during computation."
  [x! y!]
  (let [n (dim x!)
        x-mean (mean x!)
        y-mean (mean y!)
        x-centered (linear-frac! 1.0 x! (- x-mean))
        y-centered (linear-frac! 1.0 y! (- y-mean))]
    (/ (dot x-centered y-centered) (dec n)))) 

; For applications requiring preservation of input data, we provide a non-destructive version:

(defn covariance
  "Computes covariance between two vectors."
  [x y]
  (with-release [x-copy (copy x)
                 y-copy (copy y)]
    (covariance! x-copy y-copy)))

; The *covariance matrix* extends this concept to multiple variables. For a data matrix $X$ with $n$ observations and $p$ variables, the covariance matrix $\Sigma$ has elements:

; $$\Sigma_{ij} = \frac{1}{n-1}\sum_{k=1}^n (x_{ki} - \bar{x}_i)(x_{kj} - \bar{x}_j)$$

(defn matrix-covariance
  "Computes the covariance matrix for a data matrix X.
   Each row represents an observation, each column a variable."
  [X]
  (let [n (mrows X)
        p (ncols X)
        result (dge p p)
        X-mean (copy X)]
    ;; Center the data
    (loop [j 0]
      (when (< j p)
        (let [col-j (col X-mean j)
              col-mean (mean col-j)]
          (linear-frac! 1.0 col-j (- col-mean))
          (recur (inc j)))))
    ;; Compute covariances
    (loop [i 0]
      (when (< i p)
        (loop [j i]
          (when (< j p)
            (let [cov (/ (dot (col X-mean i) (col X-mean j))
                         (dec n))]
              (entry! result i j cov)
              (when (not= i j)
                (entry! result j i cov)))
            (recur (inc j))))
        (recur (inc i))))
    result))

(matrix-covariance (dge 3 2 [4 2 3 4 50 6]))

; ## Correlation Analysis

; *Correlation* normalizes *covariance* to produce a standardized measure of *linear association* between variables. The *Pearson correlation coefficient* is defined as:

; $$\rho_{xy} = \frac{cov(x,y)}{\sigma_x\sigma_y}$$

; Our implementation centers the data and normalizes by the vector *norms*:

(defn correlation!
  "Computes correlation between two vectors in-place.
   Modifies both input vectors during computation."
  [x! y!]
  (let [x-mean (mean x!)
        y-mean (mean y!)
        x-centered (linear-frac! 1.0 x! (- x-mean))
        y-centered (linear-frac! 1.0 y! (- y-mean))
        x-norm (nrm2 x-centered)
        y-norm (nrm2 y-centered)]
    (if (and (not (zero? x-norm))
             (not (zero? y-norm)))
      (/ (dot x-centered y-centered)
         (* x-norm y-norm))
      0.0)))

; As with covariance, we provide a non-destructive version:

(defn correlation
  "Computes correlation between two vectors."
  [x y]
  (with-release [x-copy (copy x)
                 y-copy (copy y)]
    (correlation! x-copy y-copy)))

; The *correlation matrix* contains all pairwise correlations between variables. For a data matrix $X$, each element is:

; $$R_{ij} = \frac{\sum_{k=1}^n (x_{ki} - \bar{x}_i)(x_{kj} - \bar{x}_j)}{\sqrt{\sum_{k=1}^n (x_{ki} - \bar{x}_i)^2\sum_{k=1}^n (x_{kj} - \bar{x}_j)^2}}$$

(defn matrix-correlation
  "Computes the correlation matrix for a data matrix X.
   Each row represents an observation, each column a variable."
  [X]
  (let [n (mrows X)
        p (ncols X)
        result (dge p p)
        X-mean (copy X)]
    ;; Center and normalize the data
    (loop [j 0]
      (when (< j p)
        (let [col-j (col X-mean j)
              col-mean (mean col-j)]
          ;; Center the column
          (linear-frac! 1.0 col-j (- col-mean))
          ;; Scale to unit norm
          (let [norm (nrm2 col-j)]
            (when (> norm 0.0)
              (scal! (/ 1.0 norm) col-j)))
          (recur (inc j)))))
    ;; Compute correlations
    (loop [i 0]
      (when (< i p)
        (loop [j i]
          (when (< j p)
            (let [corr (dot (col X-mean i) (col X-mean j))]
              (entry! result i j corr)
              (when (not= i j)
                (entry! result j i corr)))
            (recur (inc j))))
        (recur (inc i))))
    result))

; ## Data Standardization

; *Standardization* transforms variables to have zero mean and unit variance. For a vector $x$, the *z-score transformation* is:

; $$z_i = \frac{x_i - \bar{x}}{\sigma_x}$$

(defn z-score
  "Standardizes a vector to have mean 0 and standard deviation 1."
  [x]
  (let [[mu sigma] (mean-variance! (copy x))
        result (copy x)]
    (linear-frac! 1.0 result (- mu))
    (when (not (zero? sigma))
      (scal! (/ 1.0 (sqrt sigma)) result))
    result))

; For matrix data, we provide functions to center and standardize all columns:

(defn center!
  "Centers a matrix by subtracting the mean of each column."
  [a!]
  (loop [j 0]
    (when (< j (ncols a!))
      (let [col-j (col a! j)]
        (linear-frac! 1.0 col-j (- (mean col-j)))
        (recur (inc j)))))
  a!)

(defn standardize!
  "Standardizes a matrix by centering and scaling to unit variance.
   Modifies the input matrix in place.
   Returns the modified matrix."
  [a!]
  (let [n (mrows a!)
        p (ncols a!)]
    ;; First center the data
    (center! a!)
    ;; Then scale to unit variance
    (loop [j 0]
      (when (< j p)
        (let [col-j (col a! j)
              ss (dot col-j col-j)  ; sum of squares of centered data
              std-dev (sqrt (/ ss (dec n)))]
          (when (> std-dev 0.0)
            (scal! (/ 1.0 std-dev) col-j))
          (recur (inc j)))))
    a!))

; ## Data Normalization

; In data analysis and machine learning, the scale of features can significantly impact algorithm performance. When different features have vastly different ranges, they can dominate or diminish the influence of other features inappropriately. *Normalization* addresses this by transforming features to comparable scales.

; ### Min-Max Scaling

; The simplest form of *normalization* maps values to the interval $[0,1]$. For a feature vector $x$, the transformation is:

; $$x_{normalized} = \frac{x - min(x)}{max(x) - min(x)}$$

; This *linear scaling* preserves zero differences and the relative ordering of values while bounding them within $[0,1]$.
; When applied to matrices, the transformation is performed column-wise, treating each feature independently.

(defn min-max!
  "Normalizes vectors or matrices to the range [0, 1] in-place.
   For matrices, normalizes each column independently.
   Modifies the input in place.
   Returns the modified input.
   Returns columns unchanged if max equals min."
  [x!]
  (if (vector? x!)
    ; Vector case
    (let [n (dim x!)
          [min-x max-x] (loop [i 1
                               min-val (entry x! 0)
                               max-val (entry x! 0)]
                         (if (< i n)
                           (let [val (entry x! i)]
                             (recur (inc i)
                                    (min min-val val)
                                    (max max-val val)))
                           [min-val max-val]))
          range (- max-x min-x)]
      (if (zero? range)
        x!  ; Return unchanged if all values are the same
        (linear-frac! (/ 1.0 range) x! (/ (- min-x) range))))
    ; Matrix case - normalize each column
    (let [m (mrows x!)
          n (ncols x!)]
      (dotimes [j n]
        (let [col-j (col x! j)
              [min-x max-x] (loop [i 1
                                   min-val (entry col-j 0)
                                   max-val (entry col-j 0)]
                             (if (< i m)
                               (let [val (entry col-j i)]
                                 (recur (inc i)
                                        (min min-val val)
                                        (max max-val val)))
                               [min-val max-val]))
              range (- max-x min-x)]
          (when-not (zero? range)
            (linear-frac! (/ 1.0 range) col-j (/ (- min-x) range)))))
      x!)))


; ### Feature Scaling

; For algorithms sensitive to the sign of inputs or requiring symmetric ranges, scaling to [-1,1] is often more appropriate. The transformation extends min-max scaling:

; $$x_{normalized} = 2 \cdot \frac{x - min(x)}{max(x) - min(x)} - 1$$

; This centers the data around zero while maintaining relative distances between points.
; For matrices, each column is scaled independently to preserve feature-specific ranges.

(defn feature-scale!
  "Normalizes vectors or matrices to the range [-1, 1] in-place.
   For matrices, normalizes each column independently.
   Modifies the input in place.
   Returns the modified input.
   Returns columns unchanged if max equals min."
  [x!]
  (if (vector? x!)
    ; Vector case
    (let [n (dim x!)
          [min-x max-x] (loop [i 1
                               min-val (entry x! 0)
                               max-val (entry x! 0)]
                         (if (< i n)
                           (let [val (entry x! i)]
                             (recur (inc i)
                                    (min min-val val)
                                    (max max-val val)))
                           [min-val max-val]))
          range (- max-x min-x)]
      (if (zero? range)
        x!  ; Return unchanged if all values are the same
        (linear-frac! (/ 2.0 range) x! (- (/ (+ min-x max-x) range)))))
    ; Matrix case - normalize each column
    (let [m (mrows x!)
          n (ncols x!)]
      (dotimes [j n]
        (let [col-j (col x! j)
              [min-x max-x] (loop [i 1
                                   min-val (entry col-j 0)
                                   max-val (entry col-j 0)]
                             (if (< i m)
                               (let [val (entry col-j i)]
                                 (recur (inc i)
                                        (min min-val val)
                                        (max max-val val)))
                               [min-val max-val]))
              range (- max-x min-x)]
          (when-not (zero? range)
            (linear-frac! (/ 2.0 range) col-j (- (/ (+ min-x max-x) range))))))
      x!)))

; ## Conclusion

; Descriptive statistics provide the foundation for understanding and preparing data for advanced analysis. The implementations in this chapter demonstrate how to efficiently compute and transform data while adhering to numerical computing best practices.

; These tools provide essential data preprocessing capabilities for statistical analysis and machine learning applications. While what we have implemented is not exhaustive, it provides a foundation for more complex statistical computations. Next steps would be to implement robust statistics and more advanced data transformations.

; ### Robust Statistics

; While *mean* and *variance* are fundamental statistical measures, they can be sensitive to *outliers* and extreme values. *Robust statistics* provide alternative measures that are less affected by anomalous data points.

; #### Median and Quantiles

; The *median* is the middle value when data is ordered. For an odd number of observations, it is the middle value; for an even number, it is the average of the two middle values. Unlike the mean, the median is not influenced by extreme values.

; More generally, *quantiles* divide ordered data into equal-sized groups. The *p-th quantile* is the value below which a proportion $p$ of the observations fall. Common quantiles include:

; - *Quartiles* (dividing data into four parts)
; - *Deciles* (ten parts)
; - *Percentiles* (hundred parts)

; #### Robust Scale Estimates

; The *Median Absolute Deviation* (MAD) is a robust measure of variability:

; $$MAD = median(|x_i - median(x)|)$$

; The MAD is particularly useful when data contains outliers that would distort the standard deviation.

; ## Implementation Notes

; ### Performance Considerations

; For large datasets, consider:

; 1. Using *in-place operations* when input preservation isn't required
; 2. *Batching operations* to minimize memory allocation
; 3. Taking advantage of *parallel processing* capabilities where available

; The implementation uses *BLAS* (Basic Linear Algebra Subprograms) operations where possible for optimal performance on numerical computations.
