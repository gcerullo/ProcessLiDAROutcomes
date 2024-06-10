# ProcessLiDAROutcomes
Code for determining proportional megatree coverage 

#1. CalculateMegatreePerHab.R
This code calculates for each 1 ha cell within a section of the landscape which spans a harvesting gradient, what the number cells are within the 1ha that pass a given threshold height. We use different threshold definitions of megatrees (45-70 in 5 m increments).


#2. ProcessScenarioMegatreeOutcomes.R
Fit a model based on the output of 1 to determine proportional coverage of each habitat type has megatrees. Then propogate to calculate scenario-wide outcomes. 


