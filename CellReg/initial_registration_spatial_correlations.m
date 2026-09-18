function [cell_to_index_map,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations]=initial_registration_spatial_correlations(maximal_distance,spatial_correlation_threshold,spatial_footprints,centroid_locations)
% This function performs an initial cell registration across sessions
% based on a chosen spatial correlation threshold.

% Inputs:
% 1. maximal_distance
% 2. spatial_correlation_threshold
% 3. spatial_footprints
% 4. centroid_locations

% Outputs:
% 1. cell_to_index_map - list of registered cells and their index in each session
% 2. registered_cells_spatial_correlations
% 3. non_registered_cells_spatial_correlations

number_of_sessions=size(spatial_footprints,2);

% Spatial correlations are computed from a sparse representation of each
% session (see footprints_to_sparse / sparse_corr2). The list of registered
% cells only needs to remember which (session, cell) created each entry
% instead of accumulating a copy of every footprint:
sparse_footprints=cell(1,number_of_sessions);
for n=1:number_of_sessions
    sparse_footprints{n}=footprints_to_sparse(spatial_footprints{n});
end

% initializing the registration with the cells from session #1:
initial_number_of_cells=sparse_footprints{1}.n;
cell_to_index_map=zeros(initial_number_of_cells,number_of_sessions);
cell_to_index_map(:,1)=1:initial_number_of_cells;
spatial_correlation_map=zeros(initial_number_of_cells,number_of_sessions);
registered_centroid_locations=centroid_locations{1};
registered_origin=[ones(initial_number_of_cells,1) , (1:initial_number_of_cells)']; % (session, cell) of each registered entry

% allocating space:
count=0;
registered_cells_spatial_correlations=zeros(1,number_of_sessions^2*initial_number_of_cells);
non_registered_cells_spatial_correlations=zeros(1,number_of_sessions^2*initial_number_of_cells);

assigned_count=0;
non_assigned_count=0;
disp('Registering cells:');
disp('Initializing list with the cells from session #1');
display_progress_bar('Terminating previous progress bars',true)
for n=2:number_of_sessions; % registering the rest of the sessions
    display_progress_bar(['Registering cells in session #' num2str(n) ' - '],false)
    new_centroids=centroid_locations{n};
    number_of_new_cells=sparse_footprints{n}.n;
    for k=1:number_of_new_cells % for each cell
        progress_tick(k,number_of_new_cells)
        is_assigned=0;
        centroid=repmat(new_centroids(k,:),size(registered_centroid_locations,1),1);
        distance_vec=sqrt(sum((centroid-registered_centroid_locations).^2,2));
        spatial_footprints_to_check=find(distance_vec<maximal_distance);
        if ~isempty(spatial_footprints_to_check) % finding the best candidate for each cell
            corr_vec=zeros(1,length(spatial_footprints_to_check));
            % candidates are looked up by the (session, cell) that created
            % them; empty spatial footprints get a correlation of 0:
            candidate_sessions=registered_origin(spatial_footprints_to_check,1);
            candidate_cells=registered_origin(spatial_footprints_to_check,2);
            if sparse_footprints{n}.sum1(k)~=0
                for candidate_session=unique(candidate_sessions)'
                    in_session=find(candidate_sessions==candidate_session);
                    these_cells=candidate_cells(in_session);
                    is_empty=sparse_footprints{candidate_session}.sum1(these_cells)==0;
                    if any(~is_empty)
                        corr_vec(in_session(~is_empty))=sparse_corr2(sparse_footprints{n},k,sparse_footprints{candidate_session},these_cells(~is_empty));
                    end
                end
            end
            [highest_corr,highest_corr_ind]=max(corr_vec);
            if highest_corr<spatial_correlation_threshold % no registration - new cell to list
                count=count+1;
                registered_origin(initial_number_of_cells+count,:)=[n , k];
                registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
                cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                cell_to_index_map(initial_number_of_cells+count,n)=k;
                spatial_correlation_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
            else % check if there is already a registered cell
                index=spatial_footprints_to_check(highest_corr_ind);
                if cell_to_index_map(index,n)==0 % register cells together
                    cell_to_index_map(index,n)=k;
                    spatial_correlation_map(index,n)=highest_corr;
                    assigned_count=assigned_count+1;
                    is_assigned=1;
                    registered_cells_spatial_correlations(1,assigned_count)=highest_corr;
                else % there is already a registered cell
                    if highest_corr>spatial_correlation_map(index,n) % switch between cells
                        count=count+1;
                        switch_cell=cell_to_index_map(index,n);
                        switch_centroid=(centroid_locations{n}(switch_cell,:));
                        registered_origin(initial_number_of_cells+count,:)=[n , switch_cell];
                        registered_centroid_locations(initial_number_of_cells+count,:)=switch_centroid;
                        cell_to_index_map(initial_number_of_cells+count,n)=switch_cell;
                        spatial_correlation_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        cell_to_index_map(index,n)=k;
                        spatial_correlation_map(index,n)=highest_corr;
                        assigned_count=assigned_count+1;
                        is_assigned=1;
                        registered_cells_spatial_correlations(1,assigned_count)=highest_corr;
                    else % no registration - new cell to list
                        count=count+1;
                        registered_origin(initial_number_of_cells+count,:)=[n , k];
                        registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
                        cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        cell_to_index_map(initial_number_of_cells+count,n)=k;
                        spatial_correlation_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                    end
                end
            end
            if is_assigned==0 % add this value to the non registered cells
                non_assigned_count=non_assigned_count+length(spatial_footprints_to_check);
                non_registered_cells_spatial_correlations(non_assigned_count-length(spatial_footprints_to_check)+1:non_assigned_count)=corr_vec;
            else % add this value to the registered cells
                temp_corr_vec=corr_vec;
                temp_corr_vec(highest_corr_ind)=[];
                non_assigned_count=non_assigned_count+length(spatial_footprints_to_check)-1;
                non_registered_cells_spatial_correlations(non_assigned_count-length(temp_corr_vec)+1:non_assigned_count)=temp_corr_vec;
            end
        else % no candidates - new cell to list
            count=count+1;
            registered_origin(initial_number_of_cells+count,:)=[n , k];
            registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
            cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
            cell_to_index_map(initial_number_of_cells+count,n)=k;
            spatial_correlation_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
        end
    end
    display_progress_bar(' done',false)
end

registered_cells_spatial_correlations(assigned_count+1:end)=[];
non_registered_cells_spatial_correlations(non_assigned_count+1:end)=[];
non_registered_cells_spatial_correlations(non_registered_cells_spatial_correlations<0.01)=[];


end
