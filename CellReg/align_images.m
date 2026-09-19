function [spatial_footprints_corrected,centroid_locations_corrected,...
    footprints_projections_corrected,centroid_projections_corrected,...
    maximal_cross_correlation,best_translations,overlapping_area,varargout] = ...
    align_images(spatial_footprints,centroid_locations,footprints_projections,...
    centroid_projections,overlapping_area_all_sessions,microns_per_pixel,...
    reference_session_index,alignment_type,sufficient_correlation_centroids,...
    sufficient_correlation_footprints,use_parallel_processing,varargin)

% This function recieves the spatial footprints from different sessions and
% finds the optimal alignment between their FOV's. This is the first step
% after loading the data, and it includes finding the optimal overall
% translations and rotations across sessions.

% Inputs:
% 1. spatial_footprints
% 2. centroid_locations
% 3. footprints_projections
% 4. centroid_projections
% 5. overlapping_area_all_sessions - the overlapping FOV
% 6. microns_per_pixel
% 7. reference_session_index
% 8. alignment_type - 'Translations' or 'Translations and Rotations'
% 9. sufficient_correlation_centroids % correlation between the centroid locations
% 10. sufficient_correlation_footprints % correlation between the spatial footprints
% 11. use_parallel_processing -  'true' for parallel processing
% 12. varargin
%   12{1}. maximal rotation/transformation_smoothness -  if 'Translations and
%   Rotations'/'Non-rigid' is used

% Outputs:
% 1. spatial_footprints_corrected
% 2. centroid_locations_corrected
% 3. footprints_projections_corrected
% 4. centroid_projections_corrected
% 5. maximal_cross_correlation
% 6. best_translations
% 7. overlapping_area
% 8. varargout
% 8{1}. displacement_fields - for non-rigid transormation

bad_algn_sessions = [];

rotation_step=0.5; % check rotations every xx degrees
minimal_rotation=0.3; % less than this rotation in degrees does not justify rotating the cells
typical_cell_size=10; % in micrometers - determines the radius that is used for gaussfit
normalized_typical_cell_size=typical_cell_size/microns_per_pixel;

number_of_sessions=size(spatial_footprints,2);
footprint_info = get_spatial_footprints(spatial_footprints{1});

adjusted_x_size=footprint_info.size(3);
adjusted_y_size=footprint_info.size(2);

% defining the outputs:
centroid_locations_corrected=cell(size(centroid_locations));
spatial_footprints_corrected=cell(size(spatial_footprints));
footprints_projections_corrected=cell(size(footprints_projections));
centroid_projections_corrected=cell(size(centroid_projections));

centroid_locations_corrected{reference_session_index}=centroid_locations{reference_session_index};
spatial_footprints_corrected{reference_session_index}=spatial_footprints{reference_session_index};
footprints_projections_corrected{reference_session_index}=footprints_projections{reference_session_index};
centroid_projections_corrected{reference_session_index}=centroid_projections{reference_session_index};

footprint_info = get_spatial_footprints(spatial_footprints{reference_session_index});

if ~isempty(footprint_info.write2path)
    footprints_ref = footprint_info.load_footprints;
    footprint = mat_to_sparse_cell(footprints_ref.footprints);
    spatial_footprints_corrected{reference_session_index} = [footprint_info.write2path, filesep,...
        'spatial_footprints_corrected_',num2str(reference_session_index),'.mat'];
    save([footprint_info.write2path, filesep, 'spatial_footprints_corrected_',...
        num2str(reference_session_index),'.mat'],'footprint','-v7.3')
    clear footprint footprints_ref
end
center_of_FOV(1)=footprint_info.size(2)/2;
center_of_FOV(2)=footprint_info.size(3)/2;

if strcmp(alignment_type,'Translations and Rotations') % if correcting for rotations as well
    maximal_rotation=varargin{1};
    rotation_vector=zeros(1,number_of_sessions-1);
    possible_rotations=-maximal_rotation:rotation_step:maximal_rotation;
    all_rotated_projections=cell(1,number_of_sessions);
    centroid_projections_rotated=cell(size(centroid_projections));

    all_rotated_projections{reference_session_index}=footprints_projections{reference_session_index};
    centroid_projections_rotated{reference_session_index}=centroid_projections{reference_session_index};
end

maximal_cross_correlation=zeros(1,number_of_sessions-1);
if strcmp(alignment_type,'Translations and Rotations')
    best_rotations=zeros(1,number_of_sessions);
end
best_x_translations=zeros(1,number_of_sessions);
best_y_translations=zeros(1,number_of_sessions);
registration_order=setdiff(1:number_of_sessions,reference_session_index);
overlapping_area=ones(adjusted_y_size,adjusted_x_size);
overlapping_area=overlapping_area.*overlapping_area_all_sessions(:,:,reference_session_index);

% Aligning the images and cells:
display_progress_bar('Terminating previous progress bars',true)
if strcmp(alignment_type,'Non-rigid') % Non-rigid alignment:
    transformation_smoothness=varargin{1};
    best_translations=zeros(2,number_of_sessions);
    displacement_fields=zeros(number_of_sessions,adjusted_y_size,adjusted_x_size,2);
    for n=1:number_of_sessions-1
        disp(['Performing non-rigid transformation for session #' num2str(registration_order(n)) ':'])
        reference_footprints_projections_corrected=footprints_projections{reference_session_index};
        temp_footprints_projections_corrected=footprints_projections{registration_order(n)};
        % (imregdemons is kept on the CPU: its GPU implementation gives a
        % measurably different displacement field, ~0.02 pixels)
        [displacement_field,temp_footprints_projections_non_rigid_corrected]=imregdemons(temp_footprints_projections_corrected,reference_footprints_projections_corrected,'AccumulatedFieldSmoothing',transformation_smoothness);       
        footprints_projections_corrected{registration_order(n)}=temp_footprints_projections_non_rigid_corrected;
        
        footprint_info = get_spatial_footprints(spatial_footprints{registration_order(n)});
        this_session_footprints_unaligned = footprint_info.load_footprints;
        this_session_footprints_unaligned = this_session_footprints_unaligned.footprints;
        
        % each footprint is warped inside a window around its support
        % (equivalent to imwarp on the full frame, see warp_footprint_stack):
        display_progress_bar('Aligning spatial footprints: ',false)
        this_session_footprints_aligned=warp_footprint_stack(this_session_footprints_unaligned,displacement_field,true);
        display_progress_bar(' done',false)
        clear this_session_footprints_unaligned
        if ~isempty(footprint_info.write2path)
            footprint = mat_to_sparse_cell(this_session_footprints_aligned);
            spatial_footprints_corrected{registration_order(n)} = [footprint_info.write2path, filesep,...
                'spatial_footprints_corrected_',num2str(registration_order(n)),'.mat'];
            save([footprint_info.write2path, filesep, 'spatial_footprints_corrected_',...
                num2str(registration_order(n)),'.mat'],'footprint','-v7.3')
        else
            spatial_footprints_corrected{registration_order(n)}=this_session_footprints_aligned;
        end
        
        [centroid_locations_corrected(registration_order(n))]=compute_centroid_locations(spatial_footprints_corrected(registration_order(n)),microns_per_pixel);
        [centroid_projections_corrected(registration_order(n))]=compute_centroids_projections(centroid_locations_corrected(registration_order(n)),spatial_footprints_corrected(registration_order(n)));
        full_FOV_correlation=normxcorr2_gathered(footprints_projections_corrected{reference_session_index},footprints_projections_corrected{registration_order(n)});
        maximal_cross_correlation(n)=max(max(full_FOV_correlation));
        displacement_fields(registration_order(n),:,:,:)=displacement_field;
    end
    varargout=cell(1,1);
    varargout{1}=displacement_fields;
else
    display_progress_bar('Terminating previous progress bars',true)
    for n=1:number_of_sessions-1
        overlapping_area_temp=overlapping_area_all_sessions(:,:,n);
        disp(['Aligning session #' num2str(registration_order(n)) ':'])
        if strcmp(alignment_type,'Translations and Rotations')
            % Searching the rotation that maximizes the cross-correlation of
            % the centroid projections (or of the footprints projections if
            % the centroids do not correlate sufficiently):
            disp('Checking for rotations')
            temp_correlations_vector=rotation_correlations(centroid_projections{reference_session_index},centroid_projections{registration_order(n)},possible_rotations,center_of_FOV,use_parallel_processing);
            if max(temp_correlations_vector)<sufficient_correlation_centroids
                temp_correlations_vector=rotation_correlations(footprints_projections{reference_session_index},footprints_projections{registration_order(n)},possible_rotations,center_of_FOV,use_parallel_processing);
            end
            [~,ind_best_rotation]=max(temp_correlations_vector);
            
            % finding the best rotation with a gaussian fit:
            rotation_range_to_check=5; % range in degrees to check for the gaussian fit
            normalized_rotation_range_to_check=rotation_range_to_check/rotation_step;
            rotation_range=round(normalized_rotation_range_to_check);
            
            % zero padding:
            if ind_best_rotation>rotation_range && ind_best_rotation<=length(possible_rotations)-rotation_range
                localized_max_correlation=temp_correlations_vector(ind_best_rotation-rotation_range:ind_best_rotation+rotation_range);
            elseif ind_best_rotation<=rotation_range
                zero_padding_size=rotation_range-ind_best_rotation+1;
                localized_max_correlation=[zeros(1,zero_padding_size) , temp_correlations_vector(1:ind_best_rotation+rotation_range)];
            elseif ind_best_rotation>length(possible_rotations)-rotation_range
                zero_padding_size=rotation_range-(length(possible_rotations)-ind_best_rotation);
                localized_max_correlation=[temp_correlations_vector(ind_best_rotation-rotation_range:end), zeros(1,zero_padding_size)];
            end
            normalized_localized_max_correlation=localized_max_correlation-min(localized_max_correlation); % transform to zero basline
            sigma_0=0.1*rotation_range;
            [~,best_rotation_temp]=gaussfit(-rotation_range:rotation_range,normalized_localized_max_correlation./sum(normalized_localized_max_correlation),sigma_0,0);
            best_rotation=possible_rotations(ind_best_rotation)+best_rotation_temp;
            rotation_vector(n)=best_rotation;
            if abs(best_rotation)>minimal_rotation
                rotated_projections=rotate_image_interp(footprints_projections{registration_order(n)}',-best_rotation,[0 0],center_of_FOV);
                overlapping_area_temp=rotate_image_interp(overlapping_area_all_sessions(:,:,n)',-best_rotation,[0 0],center_of_FOV)';
                all_rotated_projections{registration_order(n)}=rotated_projections';
            else
                all_rotated_projections{registration_order(n)}=footprints_projections{registration_order(n)};
            end
            footprint_info = get_spatial_footprints(spatial_footprints{registration_order(n)});
            unrotated_spatial_footprints = footprint_info.load_footprints;
            unrotated_spatial_footprints = unrotated_spatial_footprints.footprints;
            
            centroid_projections_rotated(registration_order(n))=compute_centroids_projections(centroid_locations(registration_order(n)),spatial_footprints(registration_order(n)));
            
        else
            footprint_info = get_spatial_footprints(spatial_footprints{registration_order(n)});
            unrotated_spatial_footprints = footprint_info.load_footprints;
            unrotated_spatial_footprints = unrotated_spatial_footprints.footprints;
        end
        
        % Finding translations with subpixel resolution:
        if strcmp(alignment_type,'Translations and Rotations')
            full_FOV_correlation=normxcorr2_gathered(all_rotated_projections{reference_session_index},all_rotated_projections{registration_order(n)});
            cross_corr_cent=normxcorr2_gathered(centroid_projections_rotated{reference_session_index},centroid_projections_rotated{registration_order(n)});
            if max(max(cross_corr_cent))<sufficient_correlation_centroids
                cross_corr_cent=full_FOV_correlation;
            end
        else
            full_FOV_correlation=normxcorr2_gathered(footprints_projections{reference_session_index},footprints_projections{registration_order(n)});
            cross_corr_cent=normxcorr2_gathered(centroid_projections{reference_session_index},centroid_projections{registration_order(n)});
            if max(max(cross_corr_cent))<sufficient_correlation_centroids
                cross_corr_cent=full_FOV_correlation;
            end
        end
        
        cross_corr_size=size(cross_corr_cent);
        partial_FOV_correlation=full_FOV_correlation(round(cross_corr_size(1)/2-cross_corr_size(1)/6):round(cross_corr_size(1)/2+cross_corr_size(1)/6)...
            ,round(cross_corr_size(2)/2-cross_corr_size(2)/6):round(cross_corr_size(2)/2+cross_corr_size(2)/6));
        maximal_cross_correlation(n)=max(max(partial_FOV_correlation));
        
        cross_corr_partial=cross_corr_cent(round(cross_corr_size(1)/2-cross_corr_size(1)/6):round(cross_corr_size(1)/2+cross_corr_size(1)/6)...
            ,round(cross_corr_size(2)/2-cross_corr_size(2)/6):round(cross_corr_size(2)/2+cross_corr_size(2)/6));
        [~,x_ind]=max(max(cross_corr_partial));
        if strcmp(alignment_type,'Translations and Rotations')
            best_rotations(registration_order(n))=best_rotation;
        end
        [~,y_ind]=max(cross_corr_partial(:,x_ind));
        
        % finding the best translation with a gaussian fit:
        gaussian_radius=round(1.5*normalized_typical_cell_size);
        % Edge guard: for a session that does not overlap the reference, the
        % correlation peak (x_ind,y_ind) can land within gaussian_radius of the
        % cross_corr_partial border, making the sub-pixel window indices below
        % (e.g. x_ind-gaussian_radius, x_ind-1) non-positive and aborting MATLAB
        % with "Array indices must be positive integers". That crash happens
        % before the maximal_cross_correlation resemblance gate further down can
        % reject the session gracefully. Clamp the peak into the interior so the
        % sub-pixel fit is always in-bounds; a genuinely non-resembling session
        % still fails that gate and is added to bad_algn_sessions as before. For
        % a well-aligned session the peak sits near the centre of
        % cross_corr_partial, so the clamp is a no-op.
        [cross_corr_partial_height,cross_corr_partial_width]=size(cross_corr_partial);
        x_ind=min(max(x_ind,gaussian_radius+1),cross_corr_partial_width-gaussian_radius);
        y_ind=min(max(y_ind,gaussian_radius+1),cross_corr_partial_height-gaussian_radius);
        temp_corr_x=cross_corr_partial(y_ind-1:y_ind+1,x_ind-gaussian_radius:x_ind+gaussian_radius);
        sigma_0=normalized_typical_cell_size/5;
        [~,mu_x_1]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(1,:)./sum(temp_corr_x(1,:)),sigma_0,0);
        [~,mu_x_2]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(2,:)./sum(temp_corr_x(2,:)),sigma_0,0);
        [~,mu_x_3]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(3,:)./sum(temp_corr_x(2,:)),sigma_0,0);
        sub_x=mean([mu_x_1 , mu_x_2 , mu_x_2 , mu_x_3]);
        temp_corr_y=cross_corr_partial(y_ind-gaussian_radius:y_ind+gaussian_radius,x_ind-1:x_ind+1);
        [~,mu_y_1]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,1)./sum(temp_corr_y(:,1)),sigma_0,0);
        [~,mu_y_2]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,2)./sum(temp_corr_y(:,2)),sigma_0,0);
        [~,mu_y_3]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,3)./sum(temp_corr_y(:,3)),sigma_0,0);
        sub_y=mean([mu_y_1 , mu_y_2 , mu_y_2 , mu_y_3]);
        x_ind=x_ind+round(cross_corr_size(2)/2-cross_corr_size(2)/6)-1;
        y_ind=y_ind+round(cross_corr_size(1)/2-cross_corr_size(1)/6)-1;
        if abs(sub_x)>1
            warning(['X axis sub-pixel correction was ' num2str(round(100*sub_x)/100) ' for session ' num2str(registration_order(n))])
            warndlg(['X axis sub-pixel correction was ' num2str(round(100*sub_x)/100) ' for session ' num2str(registration_order(n))])
            x_ind_sub=x_ind;
        else
            x_ind_sub=x_ind+sub_x;
        end
        if abs(sub_y)>1
            warning(['Y axis sub-pixel correction was ' num2str(round(100*sub_y)/100) ' for session ' num2str(registration_order(n))])
            warndlg(['Y axis sub-pixel correction was ' num2str(round(100*sub_y)/100) ' for session ' num2str(registration_order(n))])
            y_ind_sub=y_ind;
        else
            y_ind_sub=y_ind+sub_y;
        end
        best_x_translations(registration_order(n))=(x_ind_sub-adjusted_x_size);
        best_y_translations(registration_order(n))=(y_ind_sub-adjusted_y_size);
        
        % aligning projections and centroid locations:
        if maximal_cross_correlation(n)>median(partial_FOV_correlation(:))+sufficient_correlation_footprints
            if strcmp(alignment_type,'Translations and Rotations')
                untranslated_footprints_projections=all_rotated_projections{registration_order(n)};
                untranslated_centroid_projections=centroid_projections_rotated{registration_order(n)};
            else
                untranslated_footprints_projections=footprints_projections{registration_order(n)};
                untranslated_centroid_projections=centroid_projections{registration_order(n)};
            end
            translated_projections=translate_projections(untranslated_footprints_projections',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_centroid_projections=translate_projections(untranslated_centroid_projections',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_overlapping_area=translate_projections(overlapping_area_temp',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_overlapping_area=translated_overlapping_area';
            footprints_projections_corrected{registration_order(n)}=translated_projections';
            centroid_projections_corrected{registration_order(n)}=translated_centroid_projections';
            
            % Rotating/translating each spatial footprint:
                    
            % Each footprint is rotated about the center of the FOV (when the
            % rotation exceeds minimal_rotation) and translated, keeping only
            % the pixels within a fixed radius of its centroid. This is done
            % window-by-window for the whole session by
            % transform_footprint_stack (equivalent to the per-cell
            % rotate_spatial_footprint / translate_spatial_footprint calls).
            if strcmp(alignment_type,'Translations and Rotations') && abs(best_rotation)>minimal_rotation
                display_progress_bar('Rotating and translating spatial footprints: ',false)
                cell_rotation=-best_rotation;
            else
                if strcmp(alignment_type,'Translations and Rotations') && abs(best_rotation)<=minimal_rotation
                    disp('No rotations required') % less than the minimal rotation that justifies rotating each cell - translating cells
                end
                display_progress_bar('Translating spatial footprints: ',false)
                cell_rotation=0;
            end
            aligned_spatial_footprints=transform_footprint_stack(unrotated_spatial_footprints,centroid_locations{registration_order(n)},cell_rotation,[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size],center_of_FOV,microns_per_pixel,true);
            display_progress_bar(' done',false)
            clear unrotated_spatial_footprints
            if ~isempty(footprint_info.write2path)
                footprint = mat_to_sparse_cell(aligned_spatial_footprints);
                spatial_footprints_corrected{registration_order(n)} = [footprint_info.write2path, filesep,...
                    'spatial_footprints_corrected_',num2str(registration_order(n)),'.mat'];
                save([footprint_info.write2path, filesep, 'spatial_footprints_corrected_',...
                    num2str(registration_order(n)),'.mat'],'footprint','-v7.3')
            else
                spatial_footprints_corrected{registration_order(n)}=aligned_spatial_footprints;
            end
            centroid_locations_corrected(registration_order(n))=compute_centroid_locations(spatial_footprints_corrected(registration_order(n)),microns_per_pixel);                
           
            overlapping_area=overlapping_area.*translated_overlapping_area;
        else % if no appropriate rotations/translations were found
            if strcmp(alignment_type,'Translations and Rotations') % rotating cells
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using non-rigid transformation'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using non-rigid transformation'])
            elseif strcmp(alignment_type,'Translations') % rotating cells
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using rotations'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using rotations'])
            else
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - could not align session'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - could not align session'])
            end
            bad_algn_sessions = [bad_algn_sessions, registration_order(n)];
        end
    end
    
    varargout{2} = bad_algn_sessions;
    
    best_translations=microns_per_pixel*[best_x_translations ; best_y_translations];
    if strcmp(alignment_type,'Translations and Rotations')
        best_translations=[best_translations ; best_rotations];
    end    
end

end

function correlations=rotation_correlations(reference_projection,projection,possible_rotations,center_of_FOV,use_parallel_processing)
% Maximal normalized cross-correlation between the reference projection and
% the projection rotated by each of the possible rotations. Runs on the GPU
% when enabled (fastest, no parallel pool needed), otherwise in a parfor
% when parallel processing is requested, otherwise serially.
correlations=zeros(1,length(possible_rotations));
if cellreg_use_gpu()
    reference_projection=gpuArray(reference_projection);
    projection=gpuArray(projection);
    for r=1:length(possible_rotations)
        rotated_image=rotate_image_interp(projection,possible_rotations(r),[0 0],center_of_FOV);
        cross_corr=normxcorr2(reference_projection,rotated_image);
        correlations(r)=gather(max(max(cross_corr)));
    end
elseif use_parallel_processing
    parfor r=1:length(possible_rotations)
        rotated_image=rotate_image_interp(projection,possible_rotations(r),[0 0],center_of_FOV);
        cross_corr=normxcorr2(reference_projection,rotated_image);
        correlations(r)=max(max(cross_corr));
    end
else
    display_progress_bar('Checking for rotations: ',false)
    for r=1:length(possible_rotations)
        progress_tick(r,length(possible_rotations))
        rotated_image=rotate_image_interp(projection,possible_rotations(r),[0 0],center_of_FOV);
        cross_corr=normxcorr2(reference_projection,rotated_image);
        correlations(r)=max(max(cross_corr));
    end
    display_progress_bar(' done',false);
end
end

function cross_corr=normxcorr2_gathered(template,image)
% normxcorr2 on the GPU when enabled, returned as a regular array.
cross_corr=gather(normxcorr2(maybe_gpu(template),maybe_gpu(image)));
end
