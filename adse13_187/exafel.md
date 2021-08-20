---
author: Amedeo Perazzo
title: ExaFEL Project KPP-2 Verification
---

# Science Challenge Problem Description

The ExaFEL challenge problem is to create an automated analysis pipeline for
serial femtosecond crystallography (SFX), also known as
*nanocrystallography*. Although the traditional data analysis pipeline
quantifies the diffracted Bragg spots by summation-integration (a Bragg spot
typically covers more than one pixel), the envisioned exascale algorithm will
model each pixel on the image and thus push the overall accuracy to the
desired level.

The basic workflow is envisioned as a parameter-optimization inverse
problem. Traditional crystallographic data analysis is used to determine
approximate starting values for the structure factors and geometric
factors. The nanoBragg/CCTBX software is used to forward-simulate the
diffraction images by using a parametric model of pixel intensity with which
we can estimate the posterior probability of the model within a Bayesian
framework. We then employ an iterative first derivative-based method to
compute better parameter estimates and repeat the cycle until convergence at
the maximum posterior probability.

Rapid feedback is crucial for tuning sample concentrations to achieve a
sufficient crystal hit rate, ensuring that adequate data are collected and
steering the experiment. The availability of exascale computing resources and
a HPC workflow that can handle incremental bursts of data in the analysis will
allow one to perform data analysis on the fly, providing immediate feedback on
the quality of the experimental data while determining the 3D structure of the
molecule simultaneously.

To show the scalability of the analysis pipeline, we plan to progressively
increase the fraction of the machine used for reconstruction while keeping the
number of diffraction images distributed across multiple nodes constant. The
goal is to distribute the images over an increasing number of nodes while
reducing the overall reconstruction time up to the point where the analysis
can keep up with the data collection rates (5 kHz).


Table: Challenge problem details

------------------------------------------------------------------------------------
Functional requirement     Minimum criteria
-------------------------  ---------------------------------------------------------
Physical phenomena and     Iterative estimation of crystallographic structure
associated models          factors from diffraction images using the size,
                           shape, and intensity profile of Bragg spots.

Numerical approach,        FFT, Quasi-Newton parameter estimation (limited-memory
and associated models      Fletcher-Boyden-Goldfarb-Shanno); Bayesian
                           estimation.

Simulation details         *Data ingest*: Ability to ingest diffraction images at
                           1 TB/s.

$\phantom{continue}$       *Memory*:
                           Ability to store between
                           0.1 and 10 PB of events data in
                           memory; each calculation will require between
                           $10^7$ (one run) and $10^9$ (one experiment)
                           diffraction images, and the size of an image will
                           be $O(10\,\textrm{MB})$.

$\phantom{continue}$       *Workflow*: Ability to ingest data while the
                           calculation is ongoing, ability to delegate data
                           across multiple nodes for analysis, ability to
                           exchange/average the parameter estimates across
                           nodes, ability to off-load the most computing
                           intensive tasks (e.g., x-ray tracing modeling step
                           with nanoBragg) to GPU accelerators.

Demonstration calculation  Run SFX against at least $O(10^7)$ images. If
requirements               resources (e.g., memory) are available, then a full
                           demonstration calculation will be performed by
                           running SFX against $O(10^9)$ images.
------------------------------------------------------------------------------------

# Demonstration Calculation

## Facility Requirements

*List significant I/O, workflow, and/or third party library requirements that need facility support (Note: this is outside of standard MPI and ECP-supported software technologies)*

## Input description and Execution

*Describe your problem inputs, setup, estimate of resource allocation, and runtime settings*

## Problem Artifacts

*Describe the artifacts produced during the simulation, e.g., output files, and the mechanisms that are used to process the output to verify that capabilities have been demonstrated.*

## Verification of KPP-2 Threshold

*Give evidence that*

1. *The code utilized exascale resources*
2. *The executed problem met challenge problem minimum criteria*
