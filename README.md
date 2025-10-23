<img width="150" alt="ALLEGRO Logo" src="https://github.com/AmirUCR/allegro/assets/46543443/d173addd-24ef-4532-a8b8-a902e9a8ec36">

# Introduction

[Notice: The current build may have issues compiling. Amir is working on a one line fix (to just pip install allegro-crispr) coming late October]

ALLEGRO (_<ins>Al</ins>gorithm for a <ins>L</ins>inear program <ins>E</ins>nabling <ins>G</ins>uide <ins>R</ins>NA <ins>O</ins>ptimization_) is a synthetic biology tool leveraging Google OR-Tools integer linear programming to design the smallest possible gRNA library to fulfill user-specified constraints.

- Design a Cas9 gRNA library for thousands of species simultaneously
- Flexible library design using an ensemble of options such as tracks, multiplicity, pre- and post-clustering, guide cutting efficacy prediction, and more
- Extremely fast and computationally efficient
- Written in Python, Cython, and C++
- Published in [_Nucleic Acids Research_, Volume 53, Issue 15, 28 August 2025, gkaf783](https://doi.org/10.1093/nar/gkaf783)

<img width="1034" height="251" alt="image" src="https://github.com/user-attachments/assets/4661f43b-0c4e-462e-94ff-c86b64c3817d" />

**ALLEGRO’s workflow**. **Step (1)** Given the gene sequence or the genome of hundreds to thousands of input species, ALLEGRO extracts Cas9 target sequences. **Step (2)** ALLEGRO builds and solves an (integer) linear program involving millions of variables. **Step (3)** The optimal solution of the linear program determines the sgRNA library with minimal size that covers all targets.

# Documentation
You may find the documentation for ALLEGRO at its [GitHub Wiki](https://github.com/ucrbioinfo/allegro/wiki).

# Support
If you run into any issues or have suggestions for ALLEGRO, please report them on our GitHub Issues tracker. It's the fastest way to get support and helps us improve ALLEGRO for everyone. You may also email the authors at their provided e-mail addresses on the publication directly.

# About
ALLEGRO has been developed and is maintained by <ins>Amir</ins>sadra Mohseni, and Stefano Lonardi at the University of California, Riverside.

@article{mohseni2025allegro,
    title = {Kingdom-wide CRISPR guide design with ALLEGRO},
    author = {Mohseni, Amirsadra and Nia, Reyhane Ghorbani and Tafrishi, Aida and López, Mario León and Liu, Xin-Zhan and Stajich, Jason E and Lonardi, Stefano and Wheeldon, Ian},
    journal = {Nucleic Acids Research},
    volume = {53},
    number = {15},
    pages = {gkaf783},
    year = {2025},
    doi = {10.1093/nar/gkaf783},
    url = {https://doi.org/10.1093/nar/gkaf783},
    eprint = {https://academic.oup.com/nar/article-pdf/53/15/gkaf783/64082253/gkaf783.pdf},
    publisher={Oxford University Press}
}

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15597071.svg)](https://doi.org/10.5281/zenodo.15597071)






