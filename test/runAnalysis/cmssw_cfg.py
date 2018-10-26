#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =[#'file:./pat.root'
#'/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root',
#'/store/mc/RunIIFall17MiniAODv2/ZToEE_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/8A226AAF-AC43-E811-AEF0-0CC47A4D764C.root',
#'/store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/50000/000DCB8B-2ADD-E711-9100-008CFAF35AC0.root',
#'/store/mc/PhaseIFall16MiniAOD/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/FlatPU28to62HcalNZSRAW_PhaseIFall16_exo52_90X_upgrade2017_realistic_v6_C1-v1/120000/304E419F-CC13-E711-93E9-FA163E0231A1.root',
'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/12BD4CC4-0751-E811-BCA9-0090FAA58D84.root',
]                          
process.maxEvents.input = -1
isMC = True
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold, trigger_filters, trigger_path_names, prescaled_trigger_filters, prescaled_trigger_path_names, prescaled_trigger_match_2018, trigger_match_2018

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)
#!/usr/bin/env python
Electrons = False

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
####################################
####################################
####################################

HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor #####

####################################
####################################
####################################


####################################
####################################
####################################

ZSkim = False #### Set TRUE to skim dy50to120 with a Z pt < 100 GeV #####

####################################
####################################
####################################

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2017_cff as OurSelection2017






# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
	'Our2016'  : OurSelection2016,
	'Our2017'  : OurSelection2017,
	#'OurNoIso' : OurSelectionDec2012,
	'Simple'   : OurSelection2017, # The selection cuts in the module listed here are ignored below.
	#'OurMuPrescaledNew'  : OurSelectionNew,
	#'OurMuPrescaled2012' : OurSelectionDec2012
	}

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
	# Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present.
    path_list.append(process.egmGsfElectronIDSequence)
	    
   
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut

    leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
    if len(trigger_filters)>0 and (cut_name=='Our2017' or cut_name=='Simple'):
    	leptons.trigger_filters = trigger_filters
	leptons.trigger_path_names = trigger_path_names
	leptons.prescaled_trigger_filters = prescaled_trigger_filters
	leptons.prescaled_trigger_path_names = prescaled_trigger_path_names

#    if isMC:
#	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')
    if  Electrons:
	    if cut_name == 'EmuVeto':
		    leptons.electron_muon_veto_dR = 0.1

    # Keep using old TuneP for past selections
#     if 'Dec2012' not in Selection.__file__:
#         leptons.muon_track_for_momentum = cms.string('TunePNew')
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dimuons.clone(src = cms.InputTag(allname))
	if len(trigger_filters) >  0 and (cut_name=='Our2017' or cut_name=='Simple'):
		alldil.tight_cut = trigger_match_2018
        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            #alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'MuPrescaled' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
            assert alldil.tight_cut == trigger_match
            if len(prescaled_trigger_filters)>0:
                alldil.tight_cut = prescaled_trigger_match_2018
            else:
                alldil.tight_cut = prescaled_trigger_match
     # Histos now just needs to know which leptons and dileptons to use.
      
	histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)


    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
    pobj = process.EventCounter * process.dileptonPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)



    if 'VBTF' not in cut_name and cut_name != 'Simple':
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
	for dataFilter in goodDataFiltersMiniAOD:
		#setattr(process,dataFilter 
		pobj = dataFilter * pobj


    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
    	pobj = process.PrescaleToCommonMiniAOD * pobj ####### Now it seams that there are no prescaled path ########
    path = cms.Path(pobj)
    setattr(process, pathname, path)




process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
				   dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
					met_src = cms.InputTag("slimmedMETs"),
					jet_src = cms.InputTag("slimmedJets"),
				   beamspot_src = cms.InputTag('offlineBeamSpot'),
				   vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
# 								TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
							TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
				   genEventInfo = cms.untracked.InputTag('generator'),
				   metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter"))
				   )

if isMC:
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
	obj = process.prunedMCLeptons
	obj.src = cms.InputTag('prunedGenParticles')

	from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
	process.SimpleNtupler.hardInteraction = hardInteraction
	if hasattr(process, 'pathSimple'):
		process.pathSimple *=obj * process.SimpleNtupler 
else:
	if hasattr(process, 'pathSimple'):
		process.pathSimple *= process.SimpleNtupler 
if isMC:
	switch_reco_process_name(process, "PAT") # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
print process.dumpPython()
