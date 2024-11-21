function PlotAllTracks(confined_tracks, directed_tracks, diffusion_tracks, subdiffusion_tracks)
    figure;
    hold on;

    % Plot each type of track with a different color
    PlotTracks(confined_tracks, 'green',  'Confined');
    PlotTracks(directed_tracks, 'red',  'Directed');
    PlotTracks(diffusion_tracks, 'black', 'Diffusion');
    PlotTracks(subdiffusion_tracks, 'blue',  'Subdiffusion');

    xlabel('X Position');
    ylabel('Y Position');
    title('All Trajectory Types');
    % legend('Location', 'bestoutside');
    axis tight;
    hold off;
end

function PlotTracks(tracks, color, label)
    for i = 1:length(tracks)
        track = tracks{i};
        plot(track.POSITION_X, track.POSITION_Y, color, 'LineWidth', 1);
    set(gca, 'Color', 'k');
    end
    % Plot an invisible line just for the legend
    plot(NaN, NaN, color, 'LineWidth', 2);
end