function [NCC_matrix, sorted_index] = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    %% Input parser
    p = inputParser;

    addRequired(p, 'I1');
    addRequired(p, 'I2');
    addRequired(p, 'Ftp1');
    addRequired(p, 'Ftp2');

    default_window_length = 25;
    checknum = @(x) isnumeric(x) && (mod(x, 2) ~= 0) && x>1;
    addOptional(p, 'window_length', default_window_length, checknum);

    default_min_corr = 0.95;
    addOptional(p, 'min_corr', default_min_corr, @(x) isnumeric(x) && x>=0 && x<=1);

    default_do_plot = false;
    addOptional(p, 'do_plot', default_do_plot, @(x) islogical(x));

    parse(p, I1, I2, Ftp1, Ftp2, varargin{:}); %varargin别漏了{:}, 'I1'和I1一样？
    window_length = p.Results.window_length;
    min_corr = p.Results.min_corr;
    do_plot = p.Results.do_plot;
    Im1 = double(I1);
    Im2 = double(I2);

    %% Feature preparation
    [row1, col1] = size(I1);
    [row2, col2] = size(I2);
    cut_size = ceil(window_length/2); %扣除windows_length在四周角上框不到的地方

    Ftp1(1, Ftp1(1,:)<cut_size) = 0;
    Ftp1(1, Ftp1(1,:)>(col1-cut_size)) = 0;
    Ftp1(2, Ftp1(2,:)<cut_size) = 0;
    Ftp1(2, Ftp1(2,:)>(row1-cut_size)) = 0;
    Ftp1(:, any(Ftp1==0, 1)) = [];
    
    Ftp2(1, Ftp2(1,:)<cut_size) = 0;
    Ftp2(1, Ftp2(1,:)>(col2-cut_size)) = 0;
    Ftp2(2, Ftp2(2,:)<cut_size) = 0;
    Ftp2(2, Ftp2(2,:)>(row2-cut_size)) = 0;
    Ftp2(:, any(Ftp2==0, 1)) = [];

    no_pts1 = size(Ftp1, 2);
    no_pts2 = size(Ftp2, 2);

    %% Normalization
    
    dist = floor(window_length/2);
    Mat_feat_1 = zeros(window_length^2, no_pts1);
    Mat_feat_2 = zeros(window_length^2, no_pts2);
    for i = 1:no_pts1
        norm_windows = double(I1((Ftp1(2,i)-dist):(Ftp1(2,i)+dist), (Ftp1(1,i)-dist):(Ftp1(1,i)+dist))); 
        %注意imread读入的时候图片的x与y是相反的
        mean1 = mean(norm_windows(:));
        std1 = std(norm_windows(:));
        norm_windows = (norm_windows - mean1)/std1;
        Mat_feat_1(:,i) = reshape(norm_windows, window_length^2, 1);
    end

    for i = 1:no_pts2
        norm_windows = double(I2((Ftp2(2,i)-dist):(Ftp2(2,i)+dist), (Ftp2(1,i)-dist):(Ftp2(1,i)+dist)));
        mean2 = mean(norm_windows(:));
        std2 = std(norm_windows(:));
        norm_windows = (norm_windows - mean2)/std2;
        Mat_feat_2(:,i) = reshape(norm_windows, window_length^2, 1);
    end

    %% NCC calculations
    NCC_matrix = (Mat_feat_2'*Mat_feat_1)/(window_length^2-1);
    NCC_matrix(NCC_matrix<min_corr) = 0;
    pre_sorted = reshape(NCC_matrix, no_pts2*no_pts1, 1); % size(Mat_feat_1) = (window_length^2, no_pts2)
    pre_sorted(pre_sorted == 0) = [];
    [sorted, sorted_index] = sort(pre_sorted(:), 2, 'descend'); %全部筛选别忘了: !!!

end