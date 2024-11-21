function [confined_tracks, directed_tracks, diffusion_tracks, subdiffusion_tracks, all_tracks] = ExtractAndAssignTracks(folderIdx,confined_table, directed_table, norm_table, subdif_table, XY_table, all_tracks, testedconstruct)
    
% Initialize cell arrays to store tracks for each type
    confined_tracks = ExtractTracks(confined_table, XY_table);
    directed_tracks = ExtractTracks(directed_table, XY_table);
    diffusion_tracks = ExtractTracks(norm_table, XY_table);
    subdiffusion_tracks = ExtractTracks(subdif_table, XY_table);

 % Create a structure to consolidate all track data
   
    % all_tracks.confined_tracks = confined_tracks;
    % all_tracks.directed_tracks = directed_tracks;
    % all_tracks.diffusion_tracks = diffusion_tracks;
    % all_tracks.subdiffusion_tracks = subdiffusion_tracks;
        % all_tracks.confined_tracks = [all_tracks.confined_tracks; confined_tracks];
        % all_tracks.directed_tracks = [all_tracks.directed_tracks; directed_tracks];
        % all_tracks.diffusion_tracks = [all_tracks.diffusion_tracks; diffusion_tracks];
        % all_tracks.subdiffusion_tracks = [all_tracks.subdiffusion_tracks; subdiffusion_tracks];

        
        % all_tracks.confined_tracks(end+1:end+length(confined_tracks)) = confined_tracks;
        % all_tracks.directed_tracks(end+1:end+length(directed_tracks)) = directed_tracks;
        % all_tracks.diffusion_tracks(end+1:end+length(diffusion_tracks)) = diffusion_tracks;
        % all_tracks.subdiffusion_tracks(end+1:end+length(subdiffusion_tracks)) = subdiffusion_tracks;
   
        
 all_tracks.confined_tracks.(['cell_' num2str(folderIdx),'_' (testedconstruct)]) = confined_tracks;
 all_tracks.directed_tracks.(['cell_' num2str(folderIdx),'_' (testedconstruct)]) = directed_tracks;
all_tracks.diffusion_tracks.(['cell_' num2str(folderIdx),'_' (testedconstruct)]) = diffusion_tracks;
all_tracks.subdiffusion_tracks.(['cell_' num2str(folderIdx),'_' (testedconstruct)]) = subdiffusion_tracks;

end







function tracks = ExtractTracks(type_table, XY_table, Parents_table)
    num_tracks = height(type_table);
    tracks = cell(num_tracks, 1);


    for i = 1:num_tracks
        
        track_id = type_table.ID(i);
        parentvalue = type_table.PARENT_ID(i);
        
        length_track = type_table.LENGTH(i);
        start_frame = type_table.START(i);
        end_frame = type_table.END(i);
          
               
        % Extract trajectory data for this track
        track_data = XY_table(XY_table.TrackID == parentvalue & ...
                                        XY_table.FRAME >= start_frame & ...
                                        XY_table.FRAME <= end_frame, :);
            
        % Sort by frame number
        % track_data = sortrows(track_data, 'FRAME');
       
        % Store the track data
        tracks{i} = track_data;
            
        end
    
    % Remove any empty cells (in case some tracks weren't found)
    tracks = tracks(~cellfun('isempty', tracks));
end