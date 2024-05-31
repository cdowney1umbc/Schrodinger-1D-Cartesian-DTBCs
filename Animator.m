% This script creates an animation showing the values of the wave function
% over time, and saves that animation in an mp4 file in the subdirectory
% 'Renders'

% File parameters
Subdir = 'Schrodinger-1D-Cartesian-DTBCs/Renders';
filename = 'test';
filepath = strcat(Subdir, '/', filename);

% Video parameters
fps = 30;

% The number of seconds in a simulation time unit
secs_per_sim_unit = 1;

% The total number of frames in the animation
num_frames = T * secs_per_sim_unit * fps;

% The y limits
ylims = [min(-MagPsi_max, V_min), max(MagPsi_max, V_max)];

% Create a figure
figh = figure;

% Allocate space for num_frames frames
frames(num_frames+1) = struct('cdata',[],'colormap',[]);

% For each frame of the animation
for i = 0:num_frames
    % Get the time step corresponding to the current frame
    % We are i / num_frames of the way through the animation, so the time
    % step should be this fraction of N
    n = floor(i * N / num_frames);

    % Clear the figure from the previous frame
    clf

    % Graph the real and imaginary parts and squared magnitude of the wave
    % function as well as the potential function
    hold on
    % The x axis covers the region of interest
    xlim([-L, L]);
    % The y axis covers everything within the maximum value of the
    % magnitude of the wave function
    ylim(ylims);

    xlabel('x');
    ylabel('');

    % Graph the functions
    plot(x_j, RePsi_jn(:, n+1), 'b');
    plot(x_j, ImPsi_jn(:, n+1), 'r');
    plot(x_j, MagPsi_jn(:, n+1), 'm');
    plot(x_j, Valigned_jn(:, n+1), 'g');
    hold off

    % Get and store the frame
    frames(i+1) = getframe(gcf);
end

% Open a video writer with the given filepath and fps as our framerate
Writer = VideoWriter(filepath, 'MPEG-4');
Writer.FrameRate = fps;
open(Writer)
% Write our frames to it
writeVideo(Writer, frames)
% And close the writer
close(Writer)
