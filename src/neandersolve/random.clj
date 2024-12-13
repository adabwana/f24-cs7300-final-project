(ns neandersolve.random
  (:require [uncomplicate.commons.core :refer [with-release]]
            [uncomplicate.neanderthal
             [core :refer [entry! transfer raw dim mm dia scal copy]]
             [native :refer [dge dv]]
             [math :refer [sqrt log cos sin pi]]]
            [clojure.core.reducers :as r]))

;; Random Number Generation
(def ^:private multiplier 1597)
(def ^:private increment 51749)
(def ^:private modulus 244944)
(def ^:dynamic *seed* 3)

(defn set-seed!
  "Sets the random seed for reproducible number generation."
  [seed]
  (alter-var-root #'*seed* (constantly seed)))

(defn- next-random
  "Generates the next random number using the Linear Congruential Generator."
  [prev-val]
  (mod (+ (* multiplier prev-val) increment) modulus))

(defn random-sequence
  "Creates a lazy sequence of random numbers."
  [seed]
  (iterate next-random seed))

;; Distribution Implementations
(defn uniform
  "Generates a matrix with uniform random values between min and max."
  ([rows cols]
   (uniform rows cols 0.0 1.0))
  ([rows cols min max]
   (let [range (- max min)
         matrix (dge rows cols)]
     (dotimes [i (* rows cols)]
       (entry! matrix i (+ min (* range (/ (next-random (+ i *seed*)) modulus)))))
     matrix)))

(uniform 3 2)

(defn box-muller-transform
  "Implements the Box-Muller transform for generating normal random variables."
  [u1 u2]
  (let [r (sqrt (* -2 (log u1)))
        theta (* 2 pi u2)]
    [(* r (cos theta))
     (* r (sin theta))]))

(defn normal
  "Generates a matrix with normal distribution using Box-Muller transform."
  ([rows cols]
   (normal rows cols 0.0 1.0))
  ([rows cols mean std]
   (let [matrix (dge rows cols)
         size (* rows cols)
         uniform-seq (take size (random-sequence *seed*))
         normal-pairs (map box-muller-transform
                          (take-nth 2 uniform-seq)
                          (take-nth 2 (rest uniform-seq)))]
     (doseq [i (range 0 size 2)
             :let [[n1 n2] (nth normal-pairs (quot i 2))]]
       (when (< i size) (entry! matrix i (+ mean (* std n1))))
       (when (< (inc i) size) (entry! matrix (inc i) (+ mean (* std n2)))))
     matrix)))

;; Special Matrix Patterns
(defn orthogonal
  "Generates a random orthogonal matrix using QR decomposition."
  [n]
  (let [A (normal n n)
        [Q R] (qr A)]
    (mm Q (sign (dia R)))))

(defn symmetric
  "Generates a random symmetric matrix."
  [n]
  (let [A (normal n n)]
    (mm A (trans A))))

(defn positive-definite
  "Generates a random positive definite matrix."
  [n]
  (let [A (normal n n)]
    (mm (trans A) A))) 