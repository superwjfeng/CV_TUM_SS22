addpath('../Homework1/');

% Bilder laden
Image1 = imread('sceneL.png');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('sceneR.png');
IGray2 = rgb_to_gray(Image2);

% Harris-Merkmale berechnen
features1 = harris_detector_upgrade(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
features2 = harris_detector_upgrade(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);

[NCC_matrix, sorted_index] = point_correspondence(IGray1,IGray2,features1,features2,'window_length',45,'min_corr', 0.90,'do_plot',true)