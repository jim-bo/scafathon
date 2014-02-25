scafathon
=========

**Genome scaffolding comparison framework**


## Overview
Scafathon is just a collection of scripts used to evaluate three genome scaffolding tools; SILP2, OPERA and MIP. It should be easy to extend this framework to add more tools. It relies on Quast for alignment based accuracy assesment.

***Warning: This code is only intended for use by professionals. Modification of the source code is required to set proper paths to external executables. Furthermore adding additional scaffolding tools requires manual coding.***

## Usage
The tool is divided up into several sub-programs which need to be run in-order. Use the "-h" argument after each of the following sub-commands for a description of their usage.
```python
python silp.py [sub-command] -h
```
1. *align:* aligns the paired reads (required bowtie2) installed
2. *pair:* pairs two existing SAM files
3. *meta_combine:* combines two pairs of SAM files at different percentages [metagenomic simulations]
4. *prep:* prepares selected scaffold algorithm by running all preprocessing code
5. *run:* runs scaffolding tool
6. *sim_eval:* evaluates scaffold if there exists a reference AGP with locations of each contig
7. *real_eval:* uses alignment based accuracy metrics

## Installation
This is primarily a python program, it relies on several python packages:
* numpy
* networkx

Also several seperate packages are required. These include nucmer, quast and [parabio](https://github.com/jim-bo/parabio).

## Disclaimer
This is a research tool written in a research enviroment. No support is offered and bugs may be present. Only one library size is support at this time. 

