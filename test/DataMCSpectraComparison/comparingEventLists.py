#run,lumi,event,trigger,oneMuon,twoMuon,oppSign,deltaR,mass,posPt,negPt
import numpy as np
import pandas as pd
ryan_df=pd.read_csv('ryanEventPrintOut.csv', sep=',')
ian_df=pd.read_csv('ianEventPrintOut.csv', sep=',')

#first section here just compars the files collum by collum:

ryan_run = ryan_df['run']
ian_run = ian_df['run']

ryan_lumi = ryan_df['lumi']
ian_lumi = ian_df['lumi']

ryan_event = ryan_df['event']
ian_event = ian_df['event']

ryan_trigger = ryan_df['trigger']
ian_trigger = ian_df['trigger']

ryan_oneMuon = ryan_df['oneMuon']
ian_oneMuon = ian_df['oneMuon']

ryan_twoMuon = ryan_df['twoMuon']
ian_twoMuon = ian_df['twoMuon']

ryan_oppSign = ryan_df['oppSign']
ian_oppSign = ian_df['oppSign']

ryan_deltaR = ryan_df['deltaR']
ian_deltaR = ian_df['deltaR']

ryan_mass = ryan_df['mass']
ian_mass = ian_df['mass']

ryan_posPt = ryan_df['posPt']
ian_posPt = ian_df['posPt']

ryan_negPt = ryan_df['negPt']
ian_negPt = ian_df['negPt']


run_comparison = (ryan_run != ian_run)
ne_run = run_comparison[run_comparison == True]
eq_run = run_comparison[run_comparison == False]

lumi_comparison = (ryan_lumi != ian_lumi)
ne_lumi = lumi_comparison[lumi_comparison == True]
eq_lumi = lumi_comparison[lumi_comparison == False]

event_comparison = (ryan_event != ian_event)
ne_event = event_comparison[event_comparison == True]
eq_event = event_comparison[event_comparison == False]

trigger_comparison = (ryan_trigger != ian_trigger)
ne_trigger = trigger_comparison[trigger_comparison == True]
eq_trigger = trigger_comparison[trigger_comparison == False]

oneMuon_comparison = (ryan_oneMuon != ian_oneMuon)
ne_oneMuon = oneMuon_comparison[oneMuon_comparison == True]
eq_oneMuon = oneMuon_comparison[oneMuon_comparison == False]

twoMuon_comparison = (ryan_twoMuon != ian_twoMuon)
ne_twoMuon = twoMuon_comparison[twoMuon_comparison == True]
eq_twoMuon = twoMuon_comparison[twoMuon_comparison == False]

oppSign_comparison = (ryan_oppSign != ian_oppSign)
ne_oppSign = oppSign_comparison[oppSign_comparison == True]
eq_oppSign = oppSign_comparison[oppSign_comparison == False]

deltaR_comparison = (ryan_deltaR != ian_deltaR)
ne_deltaR = deltaR_comparison[deltaR_comparison == True]
eq_deltaR = deltaR_comparison[deltaR_comparison == False]

mass_comparison = (ryan_mass != ian_mass)
ne_mass = mass_comparison[mass_comparison == True]
eq_mass = mass_comparison[mass_comparison == False]

posPt_comparison = (ryan_posPt != ian_posPt)
ne_posPt = posPt_comparison[posPt_comparison == True]
eq_posPt = posPt_comparison[posPt_comparison == False]

negPt_comparison = (ryan_negPt != ian_negPt)
ne_negPt = negPt_comparison[negPt_comparison == True]
eq_negPt = negPt_comparison[negPt_comparison == False]

print "number of run not equal: ", len(ne_run)
print "number of lumi not equal: ", len(ne_lumi)
print "number of event not equal: ", len(ne_event)
print "number of trigger not equal: ", len(ne_trigger)
print "number of oneMuon not equal: ", len(ne_oneMuon)
print "number of twoMuon not equal: ", len(ne_twoMuon)
print "number of oppSign not equal: ", len(ne_oppSign)
print "number of deltaR not equal: ", len(ne_deltaR)
print "number of mass not equal: ", len(ne_mass)
print "number of posPt not equal: ", len(ne_posPt)
print "number of negPt not equal: ", len(ne_negPt)

print "\nnumber of events are different? Let's see whats going on there. Print out the first 10 that are different:\n"


indices_different_event =  ryan_event[event_comparison == True].index.values
indices_same_event =  ryan_event[event_comparison == False].index.values

for count, index in enumerate(indices_different_event):
	print "miniAOD {} nanoAOD {}".format(ryan_event.loc[index], ian_event.loc[index])
	if count > 10: break


print "\nWell, we should probably 1. figure out why and 2. see what the events are like that don't have different events:\n"

ryan_run = ryan_run[event_comparison == False]
ian_run = ian_run[event_comparison == False]

ryan_lumi = ryan_lumi[event_comparison == False]
ian_lumi = ian_lumi[event_comparison == False]

ryan_event = ryan_event[event_comparison == False]
ian_event = ian_event[event_comparison == False]

ryan_trigger = ryan_trigger[event_comparison == False]
ian_trigger = ian_trigger[event_comparison == False]

ryan_oneMuon = ryan_oneMuon[event_comparison == False]
ian_oneMuon = ian_oneMuon[event_comparison == False]

ryan_twoMuon = ryan_twoMuon[event_comparison == False]
ian_twoMuon = ian_twoMuon[event_comparison == False]

ryan_oppSign = ryan_oppSign[event_comparison == False]
ian_oppSign = ian_oppSign[event_comparison == False]

ryan_deltaR = ryan_deltaR[event_comparison == False]
ian_deltaR = ian_deltaR[event_comparison == False]

ryan_mass = ryan_mass[event_comparison == False]
ian_mass = ian_mass[event_comparison == False]

ryan_posPt = ryan_posPt[event_comparison == False]
ian_posPt = ian_posPt[event_comparison == False]

ryan_negPt = ryan_negPt[event_comparison == False]
ian_negPt = ian_negPt[event_comparison == False]


run_comparison = (ryan_run != ian_run)
ne_run = run_comparison[run_comparison == True]
eq_run = run_comparison[run_comparison == False]

lumi_comparison = (ryan_lumi != ian_lumi)
ne_lumi = lumi_comparison[lumi_comparison == True]
eq_lumi = lumi_comparison[lumi_comparison == False]

event_comparison = (ryan_event != ian_event)
ne_event = event_comparison[event_comparison == True]
eq_event = event_comparison[event_comparison == False]

trigger_comparison = (ryan_trigger != ian_trigger)
ne_trigger = trigger_comparison[trigger_comparison == True]
eq_trigger = trigger_comparison[trigger_comparison == False]

oneMuon_comparison = (ryan_oneMuon != ian_oneMuon)
ne_oneMuon = oneMuon_comparison[oneMuon_comparison == True]
eq_oneMuon = oneMuon_comparison[oneMuon_comparison == False]

twoMuon_comparison = (ryan_twoMuon != ian_twoMuon)
ne_twoMuon = twoMuon_comparison[twoMuon_comparison == True]
eq_twoMuon = twoMuon_comparison[twoMuon_comparison == False]

oppSign_comparison = (ryan_oppSign != ian_oppSign)
ne_oppSign = oppSign_comparison[oppSign_comparison == True]
eq_oppSign = oppSign_comparison[oppSign_comparison == False]

deltaR_comparison = (ryan_deltaR != ian_deltaR)
ne_deltaR = deltaR_comparison[deltaR_comparison == True]
eq_deltaR = deltaR_comparison[deltaR_comparison == False]

mass_comparison = (ryan_mass != ian_mass)
ne_mass = mass_comparison[mass_comparison == True]
eq_mass = mass_comparison[mass_comparison == False]

posPt_comparison = (ryan_posPt != ian_posPt)
ne_posPt = posPt_comparison[posPt_comparison == True]
eq_posPt = posPt_comparison[posPt_comparison == False]

negPt_comparison = (ryan_negPt != ian_negPt)
ne_negPt = negPt_comparison[negPt_comparison == True]
eq_negPt = negPt_comparison[negPt_comparison == False]

print "total events {} after removing different event numbers {}".format(len(ryan_df),len(ryan_trigger))
print "number of run not equal: ", len(ne_run)
print "number of lumi not equal: ", len(ne_lumi)
print "number of event not equal: ", len(ne_event)
print "number of trigger not equal: ", len(ne_trigger)
print "number of oneMuon not equal: ", len(ne_oneMuon)
print "number of twoMuon not equal: ", len(ne_twoMuon)
print "number of oppSign not equal: ", len(ne_oppSign)
print "number of deltaR not equal: ", len(ne_deltaR)
print "number of mass not equal: ", len(ne_mass)
print "number of posPt not equal: ", len(ne_posPt)
print "number of negPt not equal: ", len(ne_negPt)

print "\n Ok, before we saw 155 events without the trigger being equal, now we see 150. So no real change. \n"

print "\n Let's look at the events with the same event number but not matching trigger step: \n"

differentTriggerSameEvent = ryan_event[trigger_comparison==True]

for index, event in enumerate(differentTriggerSameEvent):
	print event
	if index > 10: break

print "..."

print "\n That's pretty boring, more interesting is are all of my events passes? \n"

differentTriggerSameEvent_index = differentTriggerSameEvent.index.values

print "Printing number miniAOD's trigger fails with same events and different trigger result: {} out of {}".format( len(ryan_trigger[trigger_comparison==True][ryan_trigger !=1 ]), len(ryan_trigger[trigger_comparison==True]))

print "Printing number nanoAOD's trigger fails with same events and different trigger result: {} out of {}".format( len(ian_trigger[trigger_comparison==True][ian_trigger !=1 ]), len(ian_trigger[trigger_comparison==True]))
