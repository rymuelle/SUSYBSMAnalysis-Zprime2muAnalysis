import FWCore.ParameterSet.Config as cms

# JMTBAD to drop
muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRDPtLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0))'),
    # Follow VBTF's matching criteria.
    maxDPtRel             = cms.double(1),
    maxDeltaR             = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)
)

trigger_pt_threshold = 40
offline_pt_threshold = 45
trigger_paths = ['HLT_Mu40_eta2p1_v%i' % i for i in (9,)]
trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i && abs(userFloat("TriggerMatchEta")) < 2.1' % locals()

overall_prescale = 600
prescaled_trigger_pt_threshold = 15
prescaled_offline_pt_threshold = 20
prescaled_trigger_paths = ['HLT_Mu15_eta2p1_v%i' % i for i in (3,)]
prescaled_trigger_match = trigger_match.replace('Trigger', 'prescaledTrigger')
