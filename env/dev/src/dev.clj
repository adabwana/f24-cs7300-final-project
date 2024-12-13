(ns dev
  (:require [scicloj.clay.v2.api :as clay]))


(defn build []
  (clay/make!
   {:format              [:quarto :html]
    :book                {:title "CS7300: Final Project"}
    :subdirs-to-sync     ["notebooks" "data" "images"]
    :source-path         ["src/index.clj"
                          "src/neandersolve/descriptive.clj"
                          "src/neandersolve/eigen.clj"
                          "src/neandersolve/pca.clj"
                          "notebooks/python/pca-example.ipynb"
                          "src/neandersolve/pca-analysis.clj"
                          "notebooks/conclusion.md"]
    :base-target-path    "docs"
    ;; :live-reload true
    :clean-up-target-dir true}))

(comment
  (build))