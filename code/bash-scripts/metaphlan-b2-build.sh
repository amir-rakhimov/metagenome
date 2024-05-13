#!/usr/bin/env bash
conda activate metaphlan-tools-4.0.2
bowtie2-build --threads 12 ~/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.fna ~/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103

