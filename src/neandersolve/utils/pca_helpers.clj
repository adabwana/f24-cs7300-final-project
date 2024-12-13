(ns neandersolve.utils.pca-helpers
  (:require
   [tablecloth.api :as tc]
   [uncomplicate.commons.core :refer [let-release with-release]]
   [uncomplicate.neanderthal.core :refer [col copy dia dim entry entry! gd ge
                                          mm mm! mmt mrows mv ncols raw rk!
                                          row submatrix trans view-ge]]
   [uncomplicate.neanderthal.linalg :refer [ev!]]
   [uncomplicate.neanderthal.native :refer [dge fge]]
   [uncomplicate.neanderthal.vect-math :refer [abs! sqrt!]]))

; ============================================================
; ## Matrix centering
; ============================================================
(defn center-rows! 
  "Centers rows of a matrix by subtracting the mean of each row."
  [a!]
  (with-release [ones (entry! (raw (row a! 0)) 1.0)
                 sum-rows (mv a! ones)]
    (rk! (/ -1.0 (dim ones)) sum-rows ones a!)))


(defn center! 
  "Centers a matrix by subtracting the mean of each column."
  [a!]
  (with-release [ones (entry! (raw (col a! 0)) 1.0)
                 sum-cols (mv (trans a!) ones)]
    (rk! (/ -1.0 (dim ones)) ones sum-cols a!)))


; ============================================================
; ## PCA
; ============================================================
(defn pca 
  "Performs PCA on a matrix. May want to center the matrix first."
  [a]
  (with-release [at (trans a)
                 cov (mmt at)
                 d (gd at (mrows cov))]
    (let-release [u (ge at (mrows cov) (mrows cov))]
      (ev! cov (view-ge (dia d)) u nil)
      (mm! 1.0 u (sqrt! (abs! d))))))


; ============================================================
; ## Reconstruction
; ============================================================
(defn pca-decode 
  "Decodes a matrix from its PCA components."
  [a]
  (with-release [cov (mmt a)
                 d (gd a (mrows cov))]
    (let-release [v (ge a (mrows cov) (mrows cov))]
      (ev! cov (view-ge (dia d)) nil v)
      v)))


#_(mm (submatrix (pca-decode (copy a1)) 0 2 3 1)
    (trans (submatrix (pca (copy a1)) 0 1 2 1)))


(defn reconstruct 
  "Reconstructs a matrix from its PCA components."
  [a & {:keys [components]
        :or {components (min (mrows a) (ncols a))}}]
  (let [n-rows (mrows a)
        n-cols (ncols a)
        max-components (min (mrows a) (ncols a))
        ;; Ensure we don't exceed max possible components
        components (min components max-components)]  
    (mm (submatrix (pca-decode (copy a)) 0 (- n-rows components) n-rows components)
        (trans (submatrix (pca (copy a)) 0 (- n-cols components) n-cols components)))))