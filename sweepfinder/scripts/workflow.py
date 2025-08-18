#!/bin/env python3
import gwf
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('/faststorage/project/EcoGenetics/people/jilong/local_adaptation/sweepfinder/scripts/workflow_sources/'))
from workflow_sources import *
configs = glob.glob('./configurations/*config.y*ml')

gwf = Workflow()
for config in configs:
    gwf = sweepfinder2_workflow(config_file = config, gwf =gwf)
