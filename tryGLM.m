%% try GLM to predict NR population firing rate using other regions's unit activity
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_units.batch';
regs = ["hpc","nr","pfc","v1"];
runBatch(batch_file,@ISRPopulGLM_,{regs,'sleepm'},'verbose',true,'sessions',16);