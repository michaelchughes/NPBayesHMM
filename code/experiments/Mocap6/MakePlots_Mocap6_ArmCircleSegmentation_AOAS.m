close all;
seqNames = {'14_14', '14_06', '14_20','13_29', '13_30', '13_31'};
armcircleTsteps = { [200:240], '', '', [105:140], '', [185:215]};

jobNames = {'SM+zDD+AnnealLin', 'SM+zDD', 'SM+cDD+AnnealLin', 'SM+cDD', 'Prior'};
jobIDs = [2134873 2134874 2134875 2134876 1515159];

EXPORT_DIR = '/home/mhughes/git/BPARHMM-NEW/figs/Mocap6/';
taskID = 1;

jj =0;
for jobID = jobIDs([1 end]);
    jj = jj+1;
    ss = 0;
    for seqID = [1 4 6]
        ss = ss + 1;
        plotAlignedStateSeqSegment( jobID, taskID, seqNames{seqID}, armcircleTsteps{ seqID }  );
        
        fname = fullfile(EXPORT_DIR, sprintf('ArmCircleSegmentation_seq%d_%s', seqID, jobNames{jj}) );
        export_fig( fname, '-eps');
        fprintf( 'Exporting %s\n', fname );
    end
end
