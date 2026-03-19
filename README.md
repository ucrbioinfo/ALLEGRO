<img width="150" alt="ALLEGRO Logo" src="https://github.com/AmirUCR/allegro/assets/46543443/d173addd-24ef-4532-a8b8-a902e9a8ec36">

# Introduction
ALLEGRO is a synthetic biology tool leveraging linear programming to design the smallest possible gRNA library to fulfill user-specified constraints.

- Design a Cas gRNA library for thousands of species simultaneously
- Flexible library design using an ensemble of options such as tracks, multiplicity, pre- and post-clustering, guide cutting efficacy prediction, and more
- Extremely fast and computationally efficient
- Written in Python, Cython, and C++

**Overview of the ALLEGRO workflow.** **Step (1)** Given the gene sequence or the genome of hundreds to thousands of input species, ALLEGRO extracts Cas target sequences. **Step (2)** ALLEGRO formulates and solves an (integer) linear program including millions of variables. **Step (3)** The optimal solution of the linear program determines the guide library with minimal size that covers all targets.

# Documentation
You may find the documentation for ALLEGRO at its [GitHub Wiki](https://github.com/ucrbioinfo/allegro/wiki).

# Support
If you run into any issues or have suggestions for ALLEGRO, please report them on our GitHub Issues tracker. It's the fastest way to get support and helps us improve ALLEGRO for everyone.

# About
ALLEGRO has been developed and is maintained by <ins>Amir</ins>sadra Mohseni, and Stefano Lonardi at the University of California, Riverside.
