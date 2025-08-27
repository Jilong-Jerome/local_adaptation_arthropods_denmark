#!/bin/env python3
import gwf
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('/faststorage/project/EcoGenetics/people/jilong/local_adaptation/local_adaptation_arthropods_denmark/pop_phylo/scripts/workflow_sources/'))
from workflow_sources import *
configs = glob.glob('./configurations/*config.y*ml')

gwf = Workflow()
for config in configs:
    gwf = pop_phylo_workflow(config_file = config, gwf =gwf)