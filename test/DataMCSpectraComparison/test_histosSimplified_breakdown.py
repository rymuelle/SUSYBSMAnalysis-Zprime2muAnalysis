#!/usr/bin/env python


# Set to True in order to select e-mu and dielectron final states
Electrons = True
# Flag, setting whether we work with simulation or data
isMC =  False

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD



process.source.fileNames =[#'file:./pat.root'
#'/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/824C363B-0AC8-E611-B4A5-20CF3027A580.root',
#'file:/cms/ldap_home/hyeahyun/zp/sample/0030B9D6-72C1-E611-AE49-02163E00E602.root',
#'file:zpbb_m200_50k_01_2016_miniaod.root',
#'file:./2016-100k/zpbb_m200_50k_01_2016_miniaod.root',
#'/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/1433EFBD-F0C3-E611-A761-34E6D7BDDEC1.root'
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/34776732-BDE7-E611-B25C-0025905A48EC.root'
#'/store/user/rymuelle/zPrime/zp200.root',
#'/store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/40000/A8A02EE4-E38A-E811-8F4D-0025904E236E.root'
#'/store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/40000/A248B288-008B-E811-9124-0CC47A5450F8.root'
#'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90001/5E32398F-CA99-E611-B30C-0CC47A6C183A.root'
#'root://xrootd-cms.infn.it//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90001/5E32398F-CA99-E611-B30C-0CC47A6C183A.root'
'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90001/5E32398F-CA99-E611-B30C-0CC47A6C183A.root'
#'/store/user/rymuelle/zPrime/sample/zp_official_2016_200GeV.root'
#'./2016-100k/zpbb_m200_50k_01_2016_miniaod.root'
               ]
process.maxEvents.input = -1 # Set to a reasonable number (e.g.100) when testing locally with cmsRun
# Set global tags
'''
for fileName in process.source.fileNames:
    if "Run2016H" in fileName:
        process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14'
    elif "Run2016" in fileName:
        process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'
    else:
        # https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_PdmVMCcampai_AN3
        process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
        isMC = True
    '''    
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

#offline_pt_threshold, prescaled_offline_pt_threshold = 25,25

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.

#process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
#process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
#process.PrescaleToCommon.overall_prescale = overall_prescale
#
#process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
#process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = False
HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor on DY. If used, the k factor will be applied to ALL samples submitted. #####
    
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
#import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2011EPS_cff as OurSelection2011EPS
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).


def ntuplify(process, fill_gen_info=False):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    if fill_gen_info:
        obj = process.prunedMCLeptons
        obj.src = cms.InputTag('prunedGenParticles')

    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD_noDiLep',
                                           met_src = cms.InputTag("slimmedMETs"),
                                           jet_src = cms.InputTag("slimmedJets"),
                                           mu_src = cms.InputTag("slimmedMuons"),
                                           #jet_src = cms.InputTag("slimmedJetsPuppi"),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           TriggerResults_src = cms.InputTag('TriggerResults', '', 'HLT'), #data
                                           genEventInfo = cms.untracked.InputTag('generator'),
                                           metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter"))
                                           )
 



def for_data(process,GT):
#def for_data(process):
    print "for_data"
    #process.GlobalTag.globaltag = GT #RunH              #change line 52
    ntuplify(process)


if not isMC:
    for_data(process, '80X_dataRun2_2016SeptRepro_v6')


process.p = cms.Path(process.SimpleNtupler)