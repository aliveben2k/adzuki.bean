#usage: python geovar_freq2.py VCF_FILE POPULATION_LIST OUTPUT_FILE

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources

from geovar import *

data_path = pkg_resources.resource_filename("geovar", "data/")

# File path to the VCF File
vcf_file = sys.argv[1].format(data_path) #put_your_vcf_file_here

# File path to the population panel file
population_panel = sys.argv[2].format(data_path) #put_your_population_list_file_here

# Writing out VCF to a Frequency Table
# The extension of output file name should be *.csv format
af_df = vcf_to_freq_table(vcf_file, pop_panel=population_panel, outfile=sys.argv[3].format(data_path), minor_allele=True) #put_your_output_path_here