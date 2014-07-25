Two-State-Metapopulation-Range-Model
====================================

Steps in the pipeline:
----------------------
- query database, return raw stand table 
- filter the data (drainage) 
- state definition, return states for each plot measurement 
- run the random forest (use data file 1), return expected observation for each plot at each measurement
- data formatting, return transitions
- define the models to evaluate
    model 1: multinomial regression applied on transitions
    model 2: state dynamics model
- prepare the optimizer
- estimate parameters
- query database for temporary plots 
- state definition for temporary plots
- predict the likelihood of every observations using the random forest model
- evaluate the likelihood for each observation of the temporary sample plots given the different models