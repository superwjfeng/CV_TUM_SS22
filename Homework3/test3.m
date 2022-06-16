addpath('../Homework1/');
addpath('../Homework2/');

%% Bilder laden
Image1 = imread('sceneL.png');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('sceneR.png');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
features1 = harris_detector_upgrade(IGray1,'segment_length',9,'k',0.05,'min_dist',40,'N',50,'do_plot',false);
features2 = harris_detector_upgrade(IGray2,'segment_length',9,'k',0.05,'min_dist',40,'N',50,'do_plot',false);

%% Korrespondenzschaetzung
correspondences = point_correspondence(IGray1,IGray2,features1,features2,'window_length',25,'min_corr',0.9,'do_plot',false);

%% Fundamentalmatrix
[F_x1, F_x2, F_A, F_V] = epa(correspondences);

%% Essentielle Matrix
load('K.mat');
[E_x1, E_x2, E_A, E_V] = epa(correspondences, K);