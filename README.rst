BioFrames
===============
.. image:: https://readthedocs.org/projects/python-template/badge/?version=latest
    :target: http://python-template.readthedocs.org/en/latest/
    :alt: Documentation Status

.. image:: https://travis-ci.org/VDBWRAIR/your_project.svg
    :target: https://travis-ci.org/VDBWRAIR/your_project

.. image:: https://coveralls.io/repos/VDBWRAIR/your_project/badge.svg
    :target: https://coveralls.io/r/VDBWRAIR/your_project


This is the template for WRAIR python projects that will help you quickly setup
a base project that you can easily expand upon.

Features
--------
Currently supports .sam, .vcf, .fastq and pileup formats. 
Basic support methods for viewing & plotting quality as a matrix
filter vcf files

To Do
-----
Accept fasta files
Simplified frame merges
frame diffing
ipython notebook examples

Notes
-----
Currently importing may not be working, run nosetests with PYTHONPATH=PYTHONPATH:/home/AMED/michael.panciera/projects/biopandas/bioframes/:$PWD nosetests
