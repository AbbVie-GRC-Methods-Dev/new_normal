"""
This file imports the function that uses multiprocessing to annotate formatted MAFs with ProbSomatic metrics.
data_splits.csv is used as input for the cohort.
"""

import context

import os
from new_normal.prob_somatic_helpers import get_prob_somatic_mafs

hugo_project_dir = os.path.dirname(os.path.abspath(__file__))

# calculates prob somatic features, annotates mafs with features, and saves mafs 
get_prob_somatic_mafs(hugo_project_dir, max_cores = 100)
