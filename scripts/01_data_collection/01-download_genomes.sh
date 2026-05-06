#!/bin/bash

python download_genomes.py -i isolates.txt -o all_genomes
gunzip -f all_genomes/*.gz
        
  