#!/usr/bin/python3

from __future__ import print_function

import pandas
pandas.options.display.max_rows = 5000
from time import time


import cobra

from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion)
#from cobra.flux_analysis import gapfill
#from cobra.flux_analysis import flux_variability_analysis
#from cobra.flux_analysis import gapfill

import copy
import os
from os.path import join
import math
from math import log
import sys

import string

total = len(sys.argv)
if total == 3: 
    cmdargs = str(sys.argv)
    data_dir = os.getcwd()
    print("--------------------------------------------------------------------------------------------------------------------")
    print("Reading sbml file %s from %s" % (str(sys.argv[1]),data_dir))
else:    
    print ("Give sbml file and name of biomass production rxn (without R_) as input\n Ex: %s ccbh.sbml PAO1_Biomass " % (sys.argv[0]))
    exit()

sbml_file = str(sys.argv[1])
sbml_model=cobra.io.read_sbml_model(join(data_dir, sbml_file))
rxn_biomass = str(sys.argv[2])
rxn_biomass_info = sbml_model.reactions.get_by_id(rxn_biomass)

nreactions=len(sbml_model.reactions)
nmetabolites=len(sbml_model.metabolites)
ngenes=len(sbml_model.genes)
print("--------------------------------------------------------------------------------------------------------------------")
print("Model %s has %d metabolites, %d reactions and %d genes" % (sbml_file, nmetabolites, nreactions, ngenes))
print("Biomass RXN (%s) :: %s" % (rxn_biomass_info.name, rxn_biomass_info.reaction))
print("Objective :: %s" % (str(sbml_model.objective)))
print("--------------------------------------------------------------------------------------------------------------------")
solution=sbml_model.optimize()
#sbml_model.summary()
max_biomass_flux=solution.f

max_biomass_flux=solution.f
if max_biomass_flux < 1e-9:
    print("No biomass production")
    print("--------------------------------------------------------------------------------------------------------------------")

    print("=====================================================================================================================")
    sbml_gapfill_file = sbml_file.replace('.sbml', '_gapfill.sbml')
    print("== Using gapfilling methods to determine the smallest number of reactions needed for biomass production ==") 
    print("=====================================================================================================================")
    sbml_gapfill=cobra.io.read_sbml_model(join(data_dir, sbml_gapfill_file))
    gapfill_file = "gapfill.info"
    output_gapfill=open(gapfill_file, "w")
    output_gapfill.write("-------------------------------------------------------------------------------------------------------------------- \n")
    output_gapfill.write("  Using gapfilling methods to determine the smallest number of reactions needed for biomass production\n") 
    output_gapfill.write("-------------------------------------------------------------------------------------------------------------------- \n")
    r = cobra.flux_analysis.growMatch(sbml_model, sbml_gapfill, iterations=10)
    fbiomass_max=0.0;
    fbiomass_max_biomass=0.0;
    nentries_min=100000;
    nentries_max_biomass=100000;
    for i, entries in enumerate(r):
        ccbh_model_gapfill = sbml_model.copy()
        nentries=0
        for e in entries:
            nentries += 1
            ri = sbml_gapfill.reactions.get_by_id(e.id)
            ccbh_model_gapfill.add_reaction(ri)
        ccbh_model_gapfill.optimize()
        fccbh_gapfill = ccbh_model_gapfill.optimize().f
        fba_fluxes_gapfill=ccbh_model_gapfill.solution.x_dict
        print("-------------------------------------------------------------------------------")
        print("------ Trial %d (%d reactions): Biomass production rate = %f ------" % (i + 1, nentries, fba_fluxes_gapfill[rxn_biomass]))
        for e in entries:
            ri = sbml_gapfill.reactions.get_by_id(e.id)
            print("%s(%s) :: %s (flux=%f)" % (e.id, ri.name, ri.reaction, fba_fluxes_gapfill[e.id]))
    exit(0)
else:
    fba_fluxes=sbml_model.solution.x_dict
    print("Max. biomass production rate = %f hr-1 (dup. time = %f min)" % (max_biomass_flux , 60.0*log(2.0)/max_biomass_flux))
    fluxes_file = "Smatrix_fluxes.txt"
    with open(fluxes_file, "w") as fluxes_output:
        for f in sbml_model.reactions:
            fluxes_output.write("\"R_%s\" %.8f \"%s\"\n" % (f.id , fba_fluxes[f.id] , f.reaction))
print("-------------------------------------------------------------------------------------------------------------------- \n")
exit(0)
