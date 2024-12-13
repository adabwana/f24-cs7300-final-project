(ns neandersolve.eigen
     (:require
      [uncomplicate.fluokitten.core :refer [fmap!]]
      [uncomplicate.neanderthal
       [core :refer [axpy! col copy dia entry entry! gd 
                     ge mm mrows mv nrm2 raw scal vctr 
                     view-ge view-tr]]
       [native :refer [dge]]
       [linalg :refer [ev! org qrf trf tri!]]
       [math :refer [abs]]]))

; # Eigenvalues and Eigenvectors

; *Linear transformations* fundamentally change vectors in both magnitude and direction. However, certain special vectors maintain their direction under transformation, being only scaled by a factor. These vectors reveal intrinsic properties of the transformation that are crucial for understanding its behavior.

; ## The Eigenvalue Equation

; For a square matrix $A$, if there exists a non-zero vector $\mathbf{v}$ and scalar $\lambda$ satisfying:

; $$A \phi = \lambda \phi$$

; then $\phi$ is called an *eigenvector* and $\lambda$ its corresponding *eigenvalue*. This equation tells us that when $A$ transforms $\phi$, the result points in the same (or opposite) direction as $\phi$, scaled by $\lambda$.

; ## Power Iteration Method

; The *power iteration method* provides a simple way to find the *dominant eigenvalue* and its corresponding eigenvector. Starting with a random vector $\mathbf{x}_0$, we repeatedly apply the transformation:

; $$\mathbf{x}_{k+1} = \frac{A\mathbf{x}_k}{\|A\mathbf{x}_k\|}$$

; This process converges to the eigenvector corresponding to the largest (in magnitude) eigenvalue. The convergence rate depends on the ratio between the largest and second-largest eigenvalues.

(defn power-iteration-step
  "Performs a single step of the power iteration method.
   Returns [new-lambda new-vector]"
  [A x]
  (let [y (mv A x)
        lambda (nrm2 y)
        x-next (scal (/ 1.0 lambda) y)]
    [lambda x-next]))

(defn power-iteration
  "Implements the power method for finding the dominant eigenvalue and eigenvector.
   Returns [eigenvalue eigenvector] after convergence or max iterations."
  ([A]
   (power-iteration A 1e-10 100))
  ([A tol max-iter]
   (let [n (mrows A)
         x0 (entry! (vctr A n) 1.0)]
    ;;  (println "Matrix dimensions:" (mrows A) "x" (ncols A))
    ;;  (println "Initial vector size:" (dim x0))
    ;;  (println "Initial vector:" x0)
     (loop [x x0
            iter 0
            lambda-prev 0.0]
      ;;  (println "\nIteration:" iter)
       (let [[lambda x-next] (power-iteration-step A x)]
         ;;  (println "Current lambda:" lambda)
         ;;  (println "Lambda diff:" (abs (- lambda lambda-prev)))
         (if (or (< (abs (- lambda lambda-prev)) tol)
                 (>= iter max-iter))
           [lambda x-next]
           (recur x-next (inc iter) lambda)))))))

; Consider a simple example matrix:
(def test-matrix (dge 3 3 [4 2 3
                           2 5 1
                           3 1 6]))

; The power iteration method finds its dominant eigenvalue and vector:
(power-iteration test-matrix)

; ## QR Algorithm

; While power iteration finds the dominant eigenvalue, the *QR algorithm* can compute all eigenvalues simultaneously. The method is based on the *QR decomposition*, where a matrix $A$ is factored into $A = QR$ with $Q$ orthogonal and $R$ upper triangular.

; The algorithm proceeds iteratively:

; 1. Start with $A_0 = A$
; 2. At step $k$, compute $A_k = Q_kR_k$
; 3. Form $A_{k+1} = R_kQ_k$

; As $k$ increases, $A_k$ converges to an upper triangular matrix with eigenvalues on the diagonal. The convergence rate depends on the separation between eigenvalues.

(defn qr-step
  "Performs a single step of QR iteration.
   Returns the next matrix in the sequence."
  [A]
  (let [fact (qrf A)  ; Get QR factorization
        ;;  _ (println "QR factorization created")
        Q (org fact)   ; Get orthogonal matrix Q
        ;;  _ (println "Q matrix extracted")
        ;;  _ (println "Q:" Q)
        ; The R matrix is stored in the :or field of the factorization
        R (view-tr (:or fact) {:uplo :upper})   ; Get upper triangular matrix R directly from factorization
        ;;  _ (println "R matrix extracted")
        _ (println "R:" R)
        result (mm R Q)]  ; Compute next iteration
    result))

(defn qr-algorithm
  "Implements the QR algorithm for computing all eigenvalues.
   Returns a vector of eigenvalues after convergence."
  ([A]
   (qr-algorithm A 1e-10 100))
  ([A tol max-iter]
   (let [A0 (copy A)]
    ;;  (println "\nStarting QR algorithm")
    ;;  (println "Initial matrix:" A0)
     (loop [Ak A0
            k 0]
      ;;  (println "\nIteration:" k)
       (let [Ak+1 (qr-step Ak)
             diff (nrm2 (axpy! -1 (dia Ak+1) (dia Ak)))]
         ;;  (println "Current diagonal:" (dia Ak))
         ;;  (println "Next diagonal:" (dia Ak+1))
         ;;  (println "Difference:" diff)
         (if (or (< diff tol) (>= k max-iter))
           (dia Ak+1)
           (recur Ak+1 (inc k))))))))

; Let's examine the convergence on our test matrix:
(qr-algorithm test-matrix)

; ## Inverse Iteration

; Once we have eigenvalues, we can find their corresponding eigenvectors using *inverse iteration*. This method is based on the observation that if $\lambda$ is an eigenvalue of $A$, then $(A - \lambda I)$ is singular, and its *null space* contains the corresponding eigenvector.

; The algorithm applies power iteration to $(A - \lambda I)^{-1}$:

; 1. Start with a random vector $\mathbf{x}_0$
; 2. Solve $(A - \lambda I)\mathbf{y}_{k+1} = \mathbf{x}_k$
; 3. Set $\mathbf{x}_{k+1} = \mathbf{y}_{k+1}/\|\mathbf{y}_{k+1}\|$

; This process converges to the eigenvector corresponding to the eigenvalue closest to $\lambda$. The convergence is typically rapid when $\lambda$ is close to an actual eigenvalue.

(defn inverse-iteration
  "Implements inverse iteration for finding eigenvector given eigenvalue.
   Returns the corresponding eigenvector after convergence."
  ([A lambda]
   (inverse-iteration A lambda 1e-10 100))
  ([A lambda tol max-iter]
   (let [n (mrows A)
        ;;  _ (println "Matrix dimension:" n)
         I (dge n n)
         _ (dotimes [i n] (entry! I i i 1.0))
         ;;  _ (println "Identity matrix:" I)
         scaled-I (scal (- lambda) I)
         ;;  _ (println "Scaled identity matrix (-λI):" scaled-I)
         ;;  _ (println "Original matrix A:" A)
         A-lambda-I (axpy! 1.0 (copy A) scaled-I)  ; A - λI
         ;;  _ (println "A - λI:" A-lambda-I)
         x0 (entry! (vctr A n) 1.0)]
     (loop [x x0
            iter 0]
      ;;  (println "\nIteration:" iter)
      ;;  (println "Current x:" x)
       (let [y (copy x)
             ;;  _ (println "y before solve:" y)
             ; Create fresh LU factorization for each solve
             LU (trf (copy A-lambda-I))
             ;;  _ (println "LU matrix:" LU)
             ; Solve using tri! on fresh LU factorization
             y-next (mv (tri! LU) y)
             y-norm (nrm2 y-next)
             ;;  _ (println "y after solve:" y-next)
             ;;  _ (println "y norm:" y-norm)
             x-next (scal (/ 1.0 y-norm) y-next)]
         ;;  (println "Next x:" x-next)
         (if (or (< (nrm2 (axpy! -1 x-next x)) tol)
                 (>= iter max-iter))
           x-next
           (recur x-next (inc iter))))))))

; Using our test matrix and an eigenvalue from QR algorithm:

(def lambda1 (entry (qr-algorithm test-matrix) 0))
(def v1 (inverse-iteration test-matrix lambda1))

; We can verify this is indeed an eigenpair:

(let [Av (mv test-matrix v1)
      lambdav (scal lambda1 v1)]
  (nrm2 (axpy! -1 Av lambdav)))

; The small residual *norm* confirms that $\lambda \phi - A \phi \approx 0$, validating our computed eigenpair.

; ## Complete Eigendecomposition

; While the previous methods find individual eigenvalues and eigenvectors, many applications require the complete *eigendecomposition*. A matrix $A$ can be decomposed as:

; $$A = \Phi \Lambda \Phi^{-1}$$

; where $\Lambda$ is a diagonal matrix of eigenvalues and $\Phi$ contains the corresponding eigenvectors as columns. For *symmetric matrices*, $\Phi$ is orthogonal, making the decomposition particularly useful for numerical computations.

; The implementation handles both symmetric and general matrices efficiently:

; - For symmetric matrices, we use specialized algorithms that preserve symmetry
; - For general matrices, we handle potential complex eigenvalues
; - In both cases, eigenvalues are sorted by magnitude for consistency

(defn eigendecomposition
  "Computes complete eigendecomposition of a matrix.
   Returns [eigenvalues eigenvectors] as matrices with eigenvalues sorted in descending order."
  [A]
  (let [n (mrows A)
        symmetric? (instance? uncomplicate.neanderthal.internal.cpp.structures.RealUploMatrix A)
        ;; Create matrices based on matrix type
        [eigenvals eigenvecs]
        (if symmetric?
          ;; For symmetric matrices - use matrix directly
          (let [d (gd A n)  ;; Diagonal matrix for eigenvalues
                evecs (ge A n n)]  ;; Matrix for eigenvectors
            (ev! A (view-ge (dia d)) nil evecs)  ;; Use symmetric matrix directly
            ;; Extract diagonal values directly to vector
            [(fmap! identity (dia d)) evecs])
          ;; For general matrices
          (let [evals (raw (ge A n 2))
                evecs (ge A n n)]
            (ev! (copy A) evals nil evecs)
            ;; Extract first column (real parts) directly
            [(fmap! identity (col evals 0)) evecs]))
        ;; Create result matrices for sorted values
        sorted-vals (vctr A n)
        sorted-vecs (ge A n n)
        ;; Find indices sorted by absolute eigenvalue magnitude
        perm (vec (sort-by #(- (abs (entry eigenvals %))) (range n)))]
    ;; Copy values in sorted order using efficient operations
    (dotimes [i n]
      (let [src-idx (perm i)]
        (entry! sorted-vals i (entry eigenvals src-idx))
        ;; Use axpy! for efficient column copy
        (axpy! 1.0 (col eigenvecs src-idx) (col sorted-vecs i))))
    [sorted-vals sorted-vecs]))

(eigendecomposition test-matrix)

(defn eigendecomposition!
  "In-place version of eigendecomposition that modifies the input matrices.
   Requires pre-allocated eigenvals and eigenvecs matrices of correct dimensions."
  [A eigenvals eigenvecs]
  (ev! (copy A) eigenvals nil eigenvecs))

; ## Verification and Testing

; To ensure the correctness of our eigendecomposition, we verify that each computed pair $(\lambda, \phi)$ satisfies the eigenvalue equation $A \phi = \lambda \phi$ within numerical tolerance.

; The verification process checks:

; 1. *Eigenvalue equation*: $\|A \phi - \lambda \phi\| \approx 0$
; 2. *Vector normalization*: $\|\phi\| = 1$
; 3. For symmetric matrices, *orthogonality*: $\phi_i^T \phi_j = 0$ for $i \neq j$

(defn is-eigenpair?
  "Verifies if (lambda, v) is an eigenpair of matrix A within tolerance."
  ([A lambda v]
   (is-eigenpair? A lambda v 1e-8))
  ([A lambda v tol]
   (let [Av (mv A v)
         lambdav (scal lambda v)]
     (< (nrm2 (axpy! -1 Av lambdav)) tol))))

(require '[neandersolve.descriptive :refer
           [center! standardize! min-max! feature-scale!]])

(defn test-eigenpairs
  "Tests and prints the eigenpairs of a given matrix.
   Standardizes the matrix before computing eigendecomposition."
  [A]
  (let [[eigenvals eigenvecs] (eigendecomposition A)
        n (mrows A)]
    (println "\nTesting eigenpairs:")
    (loop [i 0]
      (when (< i n)
        (let [lambda (entry eigenvals i)
              v (col eigenvecs i)
              Av (mv A v)
              lambdav (scal lambda v)
              diff (nrm2 (axpy! -1 Av lambdav))]
          (println "Eigenvalue" i ":" lambda)
          (println "Eigenvector" i ":" (seq v))
          (println "Difference |Av - λv|:" diff)
          (println "Is eigenpair?" (is-eigenpair? A lambda v))
          (println)
          (recur (inc i)))))
    [eigenvals eigenvecs]))

; We can test our implementation with different matrix preprocessing:

(test-eigenpairs test-matrix)                           ; Original matrix
(test-eigenpairs (standardize! (copy test-matrix)))     ; Standardized
(test-eigenpairs (center! (copy test-matrix)))          ; Centered
(test-eigenpairs (min-max! (copy test-matrix)))         ; Min-max scaled
(test-eigenpairs (feature-scale! (copy test-matrix)))   ; Feature scaled

; ## Matrix Powers and Similarity

; One powerful application of eigendecomposition is efficient computation of matrix powers. For a diagonalizable matrix $A = \Phi \Lambda \Phi^{-1}$, we have:

; $$A^k = (\Phi \Lambda \Phi^{-1})^k = \Phi \Lambda^k \Phi^{-1}$$

; This allows us to compute high powers of matrices without repeated multiplication, which is particularly useful in:

; - *Markov chains*: Computing steady-state distributions
; - *Dynamical systems*: Analyzing long-term behavior
; - *Network analysis*: Computing multi-step connections

; The implementation demonstrates this with a simple example:

(comment
  (let [A (dge 2 2 [-4 3 -6 5])
        evec (dge 2 2)
        eval (ev! (copy A) (raw A) nil evec)
        inv-evec (tri! (trf evec))
        d (mm inv-evec (mm A evec))]
    (fmap! (pow 9) (dia d))
    d))

#_(defn eigendecomposition1
    "Computes complete eigendecomposition of a matrix.
    Returns [eigenvalues eigenvectors] as matrices with eigenvalues sorted in descending order."
    [A]
    (let [n (mrows A)
          ;; _ (println "Matrix size:" n)
          ;; _ (println "Input matrix A:" A)
          ; Create matrices with proper format for eigenvalues (nx2) and eigenvectors (nxn)
          eigenvals (-> A copy raw (view-ge n 2))  ; Take first 2 columns
          eigenvecs (raw A)  ; Matrix to store eigenvectors
          ;; _ (println "Created eigenvals matrix:" eigenvals)
          ;; _ (println "Created eigenvecs matrix:" eigenvecs)
          ; Compute decomposition in-place
          _ (ev! (copy A) eigenvals nil eigenvecs)
          ;; _ (println "After ev! eigenvals:" eigenvals)
          ;; _ (println "After ev! eigenvecs:" eigenvecs)
          ; Sort eigenvalues and eigenvectors
          eig-pairs (map-indexed (fn [i _] [(entry eigenvals i 0) i]) (range n))
          sorted-indices (map second (sort-by (comp - first) eig-pairs))
          ; Create new matrices for sorted results
          sorted-vals (vctr A n)
          sorted-vecs (raw A)
          ; Copy values in sorted order
          _ (doseq [[new-idx old-idx] (map-indexed vector sorted-indices)]
              (entry! sorted-vals new-idx (entry eigenvals old-idx 0))
              (axpy! 1.0 (col eigenvecs old-idx) (col sorted-vecs new-idx)))]
      [sorted-vals sorted-vecs]))

