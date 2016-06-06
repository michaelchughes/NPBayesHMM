close all;
seqNames = {'14_14', '14_06', '14_20','13_29', '13_30', '13_31'};
jogTsteps = { [42:82], '', '', '', [165:185], [145:165]};

jobNames = {'SM+zDD+AnnealLin', 'SM+zDD', 'SM+cDD+AnnealLin', 'SM+cDD', 'Prior'};
jobIDs = [2134873 2134874 2134875 2134876 1515159];

EXPORT_DIR = '/home/mhughes/git/BPARHMM-NEW/figs/Mocap6/';
taskID = 1;

jj =0;
for jobID = jobIDs( [1 end] );
    jj = jj+1;
    ss = 0;
    for seqID = [1 5 6]
        ss = ss + 1;
        plotAlignedStateSeqSegment( jobID, taskID, seqNames{seqID}, jogTsteps{ seqID }  );
        
        fname = fullfile(EXPORT_DIR, sprintf('JogSegmentation_seq%d_%s', seqID, jobNames{jj}) );
        export_fig( fname, '-eps');
        fprintf( 'Exporting %s\n', fname );
    end
end
