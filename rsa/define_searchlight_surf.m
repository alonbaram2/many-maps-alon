% RSA_wrapperMona(whatDoIDo,subj,session,design_name,type,varargin) % varargin{1} is usually RSAmask for mva_run_searchlight 

                
function define_searchlight_surf(root, subj)
        
        ROI = false;
        
        % Create directories
        slDir= fullfile(root,'rsa_alon',subj,'searchlight');
        caretDir = fullfile(root,'rsa_alon',subj,'caretDirs');
        
        % Brain mask in functional space
        RSAmask_fname = fullfile(root,subj,'bestStruct','brain_mask.nii');
        fname = 'WhBr_surf_r10_v100.mat';
        Vmask = rsa.readMask(RSAmask_fname);
        
        if ~exist(slDir,'dir') % create a folder for 1st-level results
            mkdir(slDir);
        end
        if ~exist(caretDir,'dir') % create a folder for 1st-level results
            mkdir(caretDir);
        end
        
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(['Running searchlight analysis  ',subj])
        disp(['Data will be stored here: ',slDir])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        subject_dir = getenv('SUBJECTS_DIR');
        fs_dir = getenv('FREESURFER_HOME');
        
        disp(['Freesurfer directory: ',fs_dir])
        
        if ~exist(fullfile(subject_dir, ['x', subj]) ,'dir')
            
            % register the right hemisphere to the left (xhemi) and register
            % both to fsaverage_sym (atlas).
            % First inverts the right hemi using 'xhemireg'.
            % Then creates registration files using 'surfreg':
            % Left hemi: lh.fsaverage_sym.sphere.reg in subject/surf
            % Right Hemi: lh.fsaverage_sym.sphere.reg in subject/xhemi/surf
            freesurfer_registerXhem({subj},subject_dir,hemisphere',[1 2]);
            
            % resample the registered subject's surfaces into atlas space (ico7). Creates new
            % remapped surfaces (.white, .pial, .inflated) and curves (.cure, .sulc, .area)
            % which will be in $SUBJECTS_DIR/['x' subj] dir (both hemisphere). Still not quite sure if new
            % surfaces are in fsaverage_sym or ico7 (fsaverage) space.
            % Effectively this means changeing the number and order of verteces of the input
            % surfaces to match ico7, but keep the information about the
            % the RAS coordinates of the nearest neighbour vertex
            % of in the (moved) subject surface to each ico vertex.
            % the topology (surface faces) is taken from ico.
            freesurfer_mapicosahedron_xhem(subj,subject_dir,'smoothing',1,'hemisphere',[1:2]);
            
            
            % Move the surfaces vertex coords to correspond to World
            % rather than RAS_tkreg coords (256x256x256,
            % like in $SUBJECTS_DIR/$sub/orig.mgz).
            % save in Caret format.
            caret_importfreesurfer(['x', subj],subject_dir,caretDir);
        end
        
        LcaretDir = fullfile(caretDir,sprintf(['x', subj]),'LeftHem');
        RcaretDir = fullfile(caretDir,sprintf(['x', subj]),'RightHem');
        white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
        pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
        topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
        
        % S is a structure with the subjects surfaces in world
        % (mm) coords. I.e. if you save S as a surface and load it in
        % Freeview, the RAS coordinates of each vertex will correspond
        % exactly to the "scanner anatomical" coords in FSLeyes.
        S         = rsa.readSurf(white,pial,topo);
        % save(fullfile(slDir,'S.mat'),'S'); % save somewhere
        % save(fullfile(slDir,'Vmask.mat'),'Vmask'); % save somewhere
        
        if ROI
            L = rsa.defineSearchlight_surface(S,Vmask,'sphere',999);       % look in radius of 999 - too big so will be whole ROI
        else
            L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[10, 100]);
        end
        save(fullfile(slDir,fname),'-struct','L'); % will save the fields of structure L   
end