#!/bin/env python3
import gwf
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('./workflow_sources/'))
from workflow_sources import *

gwf = Workflow()
configs = glob.glob('./configurations/*config.y*ml')

for config in configs:
    gwf = lr_distribution_workflow(config_file = config, gwf = gwf)