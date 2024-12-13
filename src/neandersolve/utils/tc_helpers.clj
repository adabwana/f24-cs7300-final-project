(ns neandersolve.utils.tc-helpers 
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
; ## Matrix/Dataset conversion
; ============================================================
(defn matrix->dataset
  "Transforms a Neanderthal matrix into a Tablecloth dataset.
   Columns are named x1 to xp, where p is the number of columns in the matrix."
  [matrix]
  (let [rows (mrows matrix)
        cols (ncols matrix)
        column-names (mapv #(keyword (str "x" (inc %))) (range cols))]
    (loop [i 0
           result (transient [])]
      (if (< i rows)
        (recur 
         (inc i)
         (conj! result 
                (loop [j 0
                       row-map (transient {})]
                  (if (< j cols)
                    (recur 
                     (inc j)
                     (assoc! row-map 
                            (nth column-names j)
                            (entry matrix i j)))
                    (persistent! row-map)))))
        (tc/dataset (persistent! result))))))

(matrix->dataset 
 (dge 3 3 [1 2 3 4 5 6 7 8 9] {:layout :row}))

(matrix->dataset 
 (fge 3 3 [1 2 3 4 5 6 7 8 9] {:layout :column}))


(defn dataset->matrix
  "Converts a Tablecloth dataset to a Neanderthal matrix."
  [dataset]
  (let [data (tc/rows dataset)]
    (fge (count data) (count (first data))
         (flatten data) {:layout :row})))

(dataset->matrix 
 (tc/dataset {:x [1 2 3]
              :y [4 5 6]}))