function [] = ResumeExperiment( jobID, taskID, Niter )
%INPUT
%  jobID : integer/name of job
%  taskID : integer/name of task
%  Niter  : # iterations of MCMC


resumeBPHMM( {jobID, taskID}, {'Niter', Niter}, {'jobID', jobID, 'taskID', taskID} );
