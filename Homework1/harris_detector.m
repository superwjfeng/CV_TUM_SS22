function [segment_length, k, tau, do_plot] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.

    %% Input parser
    p = inputParser;
    addRequired(p, 'input_image');
    
    segment_length_default = 15;
    addOptional(p, 'segment_length', segment_length_default, @(x) isnumeric(x) && (x > 1) && (mod(x,2)~=0));

    

    
    
end