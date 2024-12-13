# Statistical Computing in Clojure: Functional Approaches to Unsupervised Learning

**Author:** Jaryt Salvo  
**Course:** CS 7300 Unsupervised Learning  
**Term:** Fall 2024

*************

## **[Project Documentation and Results](https://adabwana.github.io/f24-cs7300-final-project/)**

## Project Overview

This work implements fundamental *statistical algorithms* using **functional programming** paradigms. The implementation leverages `Clojure`'s *immutable data structures* and **pure functions** to develop a mathematically rigorous framework for *numerical computing*. Our approach emphasizes *computational efficiency* while maintaining *mathematical precision* through careful algorithm selection and implementation.

### Core Implementation

The project consists of three primary computational modules:

**1. Statistical Foundations**

- Implementation of **numerically stable algorithms** for *variance* and *covariance* computation
- Robust *central tendency measures* optimized for large-scale data processing
- **Matrix operations** designed for efficient *statistical analysis*

**2. Eigenvalue Decomposition**

- **Power iteration** methods for computing *dominant eigenvalues*
- Complete **eigendecomposition** through `QR` algorithm implementation
- Specialized **inverse iteration** techniques for *eigenvector* computation

**3. Principal Component Analysis**

- *Covariance*-based **PCA** implementation with mathematical rigor
- Optimized **matrix transformations** for *dimensionality reduction*
- Comparative analysis against established implementations

### Technical Implementation

The codebase utilizes modern `Clojure` libraries and practices:

- `Neanderthal` for **high-performance numerical computations**
- **Systematic validation** through *property-based testing*
- Integrated documentation with *mathematical derivations*

## Development Environment

The project uses **containerized development environments** to ensure *computational reproducibility* across systems. This approach eliminates environment-specific issues and maintains consistent behavior across different operating systems.

### System Requirements
- **Docker Desktop**
- **Visual Studio Code** with `Dev Containers` extension

### Environment Setup

1. Install the required software:
   - [Docker Desktop](https://docs.docker.com/get-docker/)
   - [Visual Studio Code](https://code.visualstudio.com/)
   - [`Dev Containers` Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

2. Clone and initialize the repository:
   ```bash
   git clone https://github.com/adabwana/f24-cs7300-final-project.git
   cd f24-cs7300-final-project
   code .
   ```

3. Select development environment:
   - Primary Container: `Clojure` environment for core implementation
   - Validation Container: `Python` environment for comparative analysis

The container configuration automatically:
- Configures language-specific toolchains
- Installs project dependencies
- Establishes consistent development environments
- Isolates project-specific configurations

This **containerized approach** ensures that all *computational results* are reproducible and that the development environment remains consistent across different systems and users.
