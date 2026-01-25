# kmed: Distance-Based k-Medoids Clustering

[![CRAN Version](https://www.r-pkg.org/badges/version/kmed)](https://CRAN.R-project.org/package=kmed)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)

## Overview

The `kmed` package provides a specialized suite of algorithms for **distance-based k-medoids clustering**. Unlike the standard k-means approach, k-medoids is highly robust against outliers and is capable of processing any distance or dissimilarity matrix. This makes it an essential tool for datasets where Euclidean distance is not applicable, such as those containing categorical, ordinal, or mixed-type variables.

The package is designed for statistical researchers and practitioners who require flexible clustering solutions for complex, non-Euclidean data structures.

## Key Features

* **Diverse Algorithms:** Implements several k-medoids variants, including:
    * **Simple K-Medoids (SKM)** and **Fast K-Medoids (FKM)**.
    * **Ranked K-Medoids (RKM)**.
    * **Increasing Number of Clusters (INCO)** for determining optimal cluster partitions.
* **Mixed-Type Data Support:** Optimized functions for calculating distances in heterogeneous datasets using specialized metrics:
    * **Gower**, **Podani**, and **Wishart** distances.
    * **Huang**, **Harikumar-PV**, and **Ahmad-Dey** dissimilarities.
* **Cluster Validation:** Built-in tools for evaluating clustering quality, including:
    * **Internal Criteria:** Silhouette Index and Shadow Values.
    * **Relative Criteria:** Bootstrap procedures for cluster stability.
* **Visualization:** Advanced support for bootstrap-based heatmaps, marked barplots, and PCA biplots to facilitate the interpretation of clustering results.

## Installation

The stable version of `kmed` is available on [CRAN](https://CRAN.R-project.org/package=kmed). You can install it directly via the R console using the following command:

```r
install.packages("kmed")
