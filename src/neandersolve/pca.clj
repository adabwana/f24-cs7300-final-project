(ns neandersolve.pca
  (:require
   [neandersolve
    [descriptive :as desc]
    [eigen :as eigen]]
   [uncomplicate.commons.core :refer [with-release]]
   [uncomplicate.fluokitten.core :refer [fmap!]]
   [uncomplicate.neanderthal
    [core :refer [axpy col copy entry entry! mm mrows mv! ncols raw rk! scal!
                  sum trans vctr]]
    [native :refer [fge]]
    [math :refer [sqrt]]
    [random :refer [rand-normal! rand-uniform!]]]))

; # Principal Component Analysis

; High-dimensional data often contains redundant or correlated features. While each feature may carry information, the true patterns often lie in lower-dimensional subspaces. *Principal Component Analysis* (PCA) provides a mathematical framework for discovering these intrinsic patterns by transforming data into a new coordinate system aligned with directions of maximum variance.

; ## Mathematical Foundation

; For a dataset with $n$ observations and $p$ features, represented as an $n \times p$ matrix $\mathbf{X}$, PCA seeks a transformation that reveals the underlying structure. This transformation comes from the *eigendecomposition* of the *covariance matrix*:

; $$\mathbf{C} = \frac{1}{n-1}\mathbf{X}^T\mathbf{X}$$

; The *eigenvectors* of $\mathbf{C}$ form the *principal components*, while their corresponding *eigenvalues* indicate the amount of variance explained by each component.

; ## Data Preprocessing

; Before computing principal components, we must ensure our data is properly *centered* and optionally *scaled*. This preprocessing step is crucial for PCA's effectiveness:

(defn center-data
  "Centers the data matrix by subtracting column means.
   Returns [centered-data column-means]."
  [X]
  (let [n (mrows X)
        means (entry! (vctr X (ncols X)) 0.0)
        centered (copy X)
        ones (entry! (vctr X n) 1.0)]
    ;; Compute means efficiently using matrix-vector multiplication
    (mv! (/ 1.0 n) (trans X) ones 0.0 means)
    ;; Center using rank-1 update with outer product of ones and means
    (rk! -1.0 ones means centered)
    [centered means]))

; The *centering operation* ensures that each feature has zero mean, removing location effects that could bias our variance calculations. For features measured on different scales, we also provide scaling to unit variance:

(defn scale-data
  "Scales the data matrix to unit variance.
   Returns [scaled-data column-stds]."
  [X]
  (let [p (ncols X)
        stds (entry! (vctr X p) 0.0)
        scaled (copy X)]
    (loop [j 0]
      (when (< j p)
        (let [col-std (sqrt (desc/variance (col X j)))]
          (entry! stds j col-std)
          (scal! (/ 1.0 col-std) (col scaled j))
          (recur (inc j)))))
    [scaled stds]))

; ## Covariance Computation

; The *covariance matrix* captures the relationships between features. For centered data $\mathbf{X}$, we compute it efficiently through matrix multiplication:

; $$\mathbf{C} = \frac{1}{n-1}\mathbf{X}^T\mathbf{X}$$

(defn compute-covariance
  "Computes the covariance matrix of centered data.
   Each row represents an observation, each column a variable."
  [X-centered]
  (let [n (mrows X-centered)
        covar (mm (trans X-centered) X-centered)]
    (scal! (/ 1.0 (dec n)) covar)))

; ## Model Fitting

; The core PCA algorithm combines preprocessing, covariance computation, and eigendecomposition into a unified workflow. This process reveals the *principal components* and their relative importance in explaining data variance:

; 1. Optionally *center* the data: $\mathbf{X}_c = \mathbf{X} - \mathbf{1}\boldsymbol{\mu}^T$
; 2. Optionally *scale* to unit variance: $\mathbf{X}_s = \mathbf{X}_c\mathbf{D}^{-1}$
; 3. Compute *covariance matrix*: $\mathbf{C} = \frac{1}{n-1}\mathbf{X}_s^T\mathbf{X}_s$
; 4. Perform *eigendecomposition*: $\mathbf{C} = \Phi\mathbf{\Lambda}\Phi^T$

(defn pca-fit
  "Fits PCA model to data matrix X.
   Returns map containing principal components, explained variance, etc."
  ([X]
   (pca-fit X true true))
  ([X center?]
   (pca-fit X center? false))
  ([X center? scale?]
   (let [[X-processed means] (if center? 
                              (center-data X)
                              [X (entry! (vctr X (ncols X)) 0.0)])
         [X-final stds] (if scale? 
                         (scale-data X-processed)
                         [X-processed (entry! (vctr X-processed (ncols X-processed)) 1.0)])
         cov-matrix (compute-covariance X-final)
         [eigenvals eigenvecs] (eigen/eigendecomposition cov-matrix)
         total-var (sum eigenvals)
         explained-var-ratio (fmap! #(/ % total-var) (copy eigenvals))]
     {:components eigenvecs
      :explained_variance eigenvals
      :explained_variance_ratio explained-var-ratio
      :means means
      :scale stds})))

; ## Data Transformation

; Once we have fitted a PCA model, we can transform new data into the *principal component space*. This transformation involves:

; 1. *Centering*: $\mathbf{X}_c = \mathbf{X} - \mathbf{1}\boldsymbol{\mu}^T$
; 2. *Scaling*: $\mathbf{X}_s = \mathbf{X}_c\mathbf{D}^{-1}$
; 3. *Projection*: $\mathbf{X}_{pca} = \mathbf{X}_s\Phi_k$

; where $\Phi_k$ contains the first $k$ principal components.

(defn transform
  "Transforms data using fitted PCA model.
   Optional n-components parameter for dimensionality reduction."
  ([X pca-model]
   (transform X pca-model (ncols (:components pca-model))))
  ([X pca-model n-components]
   (let [means (:means pca-model)
         scale (:scale pca-model)
         ;; Center the data
         X-centered (if (some? means)
                     (let [means-mat (raw X)]
                       (loop [i 0]
                         (when (< i (mrows X))
                           (loop [j 0]
                             (when (< j (ncols X))
                               (entry! means-mat i j (entry means j))
                               (recur (inc j))))
                           (recur (inc i))))
                       (axpy -1.0 means-mat X))
                     X)
         ;; Scale the centered data
         X-scaled (if (some? scale)
                   (let [scaled (copy X-centered)]
                     (loop [j 0]
                       (when (< j (ncols scaled))
                         (scal! (/ 1.0 (entry scale j)) (col scaled j))
                         (recur (inc j))))
                     scaled)
                   X-centered)
         ;; Select components for dimensionality reduction
         components (let [full-components (:components pca-model)
                          selected (raw (mm (trans X) X))]
                     (loop [i 0]
                       (when (< i (mrows full-components))
                         (loop [j 0]
                           (when (< j n-components)
                             (entry! selected i j (entry full-components i j))
                             (recur (inc j))))
                         (recur (inc i))))
                     selected)]
     ;; Transform the data using the principal components
     (mm X-scaled components))))

(transform (fge 2 3 [1 2 3 4 5 6]) (pca-fit (fge 2 3 [1 2 3 4 5 6])))

; ## Inverse Transformation

; To reconstruct data from its principal component representation, we reverse the transformation process:

; 1. *Back-projection*: $\mathbf{X}_s = \mathbf{X}_{pca}\Phi_k^T$
; 2. *Unscaling*: $\mathbf{X}_c = \mathbf{X}_s\mathbf{D}$
; 3. *Uncentering*: $\mathbf{X} = \mathbf{X}_c + \mathbf{1}\boldsymbol{\mu}^T$

(defn inverse-transform
  "Reconstructs original data from transformed data."
  [X-transformed pca-model]
  (let [;; Project back to original feature space
        components (let [full-components (:components pca-model)
                         selected (raw (mm (trans X-transformed) X-transformed))]
                     (loop [i 0]
                       (when (< i (mrows full-components))
                         (loop [j 0]
                          (when (< j (ncols X-transformed))
                            (entry! selected i j (entry full-components i j))
                            (recur (inc j))))
                         (recur (inc i))))
                     selected)
        X-scaled (mm X-transformed (trans components))
        ;; Unscale the data
        X-unscaled (if (some? (:scale pca-model))
                    (let [scaled (copy X-scaled)]
                      (loop [j 0]
                        (when (< j (ncols scaled))
                          (scal! (entry (:scale pca-model) j) (col scaled j))
                          (recur (inc j))))
                      scaled)
                    X-scaled)
        ;; Uncenter the data
        result (if (some? (:means pca-model))
                (let [means-mat (raw X-unscaled)]
                  (loop [i 0]
                    (when (< i (mrows X-unscaled))
                      (loop [j 0]
                        (when (< j (ncols X-unscaled))
                          (entry! means-mat i j (entry (:means pca-model) j))
                          (recur (inc j))))
                      (recur (inc i))))
                  (axpy 1.0 means-mat X-unscaled))
                X-unscaled)]
    result))

; ## Utility Functions

; ### Explained Variance Analysis

; The *explained variance ratio* helps determine the optimal number of components to retain. The *cumulative explained variance* shows how much of the total variance is captured by the first $k$ components:

(defn explained-variance-cumsum
  "Computes cumulative explained variance ratio."
  [pca-model]
  (reductions + (:explained_variance_ratio pca-model)))

(defn n-components-for-variance
  "Determines number of components needed to explain desired variance ratio."
  [pca-model target-variance]
  (let [cumsum (explained-variance-cumsum pca-model)]
    (inc (count (take-while #(< % target-variance) cumsum))))) 

; ### Example Usage

; Transform a simple 2x3 matrix:
(transform (fge 2 3 [1 2 3 4 5 6]) 
          (pca-fit (fge 2 3 [1 2 3 4 5 6])))

; Reconstruct the original data:
(inverse-transform (fge 2 3 [1.22 -1.22 0 0 0 0 0]) 
                  (pca-fit (fge 2 3 [1 2 3 4 5 6])))


#_(with-release [a (rand-normal! (fge 100 10000))]
    (time (pca-fit a)))

#_(defn inverse-transform
    "Reconstructs original data from transformed data."
    [X-transformed pca-model]
    (let [;; Project back to original feature space using transposed components
          components (let [full-components (:components pca-model)]
                          selected (raw (mmt X-transformed))
                      (loop [i 0]
                        (when (< i (mrows full-components))
                          (loop [j 0]
                            (when (< j (ncols X-transformed))
                              (entry! selected i j (entry full-components i j))
                              (recur (inc j))))
                          (recur (inc i))))
                      selected)
          ;; Use in-place matrix multiplication
          X-scaled (mm! 1.0 X-transformed (trans components) (raw (mmt X-transformed)))
          ;; Unscale efficiently
          X-unscaled (if (some? (:scale pca-model))
                      (let [scaled (copy X-scaled)]
                        (loop [j 0]
                          (when (< j (ncols scaled))
                            (scal! (entry (:scale pca-model) j) (col scaled j))
                            (recur (inc j))))
                        scaled)
                      X-scaled)
          ;; Uncenter using rank-1 update
          result (if (some? (:means pca-model))
                  (let [means-mat (raw X-unscaled)
                        ones (entry! (raw (row X-unscaled 0)) 1.0)]
                    (rk! (/ 1.0 (dim ones)) (:means pca-model) ones means-mat)
                    (axpy 1.0 means-mat X-unscaled))
                  X-unscaled)]
      result))
