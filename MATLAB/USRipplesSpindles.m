%% 1. ripple-spindle coupling
[mat,time_bins] = runBatch(batch_file,@RipSpinPETH_,{},verbose=true);

mat_pool = cellfun(@(x) x.ref,mat(~cellfun(@isempty, mat)),'UniformOutput',false);
mat_pool = vertcat(mat_pool{:});
mat_inv_pool = cellfun(@(x) x.inverse,mat(~cellfun(@isempty, mat)),'UniformOutput',false);
mat_inv_pool = vertcat(mat_inv_pool{:});
mat_on_pool = cellfun(@(x) x.on,mat(~cellfun(@isempty, mat)),'UniformOutput',false);
mat_on_pool = vertcat(mat_on_pool{:});
mat_off_pool = cellfun(@(x) x.off,mat(~cellfun(@isempty, mat)),'UniformOutput',false);
mat_off_pool = vertcat(mat_off_pool{:});


% CHECK THAT:
% 1. you are correctly excluding spindles with PETH (after restricting rips)
% 2. effect is maybe between US / no US

makeFigure('ref');
semplot(time_bins{2},mat_pool,myColors(3),'smooth',1.5,'patchProp',{'FaceAlpha',0.3})
xlabel('time from ripple peak'), ylabel('spindle occurence (a.u.)')

makeFigure('spin_us');
semplot(time_bins{2},mat_on_pool,myColors(2),'smooth',1.5,'patchProp',{'FaceAlpha',0.3})
semplot(time_bins{2},mat_off_pool,myColors(1),'smooth',1.5,'patchProp',{'FaceAlpha',0.3})
xlabel('time from ripple peak'), ylabel('spindle occurence (a.u.)')

makeFigure('spin')
semplot(time_bins{2},mat_inv_pool,myColors(1),'smooth',1.5,'patchProp',{'FaceAlpha',0.3})
xlabel('time from spindle peak'), ylabel('ripple occurence (a.u.)')