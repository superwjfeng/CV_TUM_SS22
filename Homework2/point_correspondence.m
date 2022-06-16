function cor = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
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

    %% NCC(Normalized Cross Correlation) calculations
%     这么写是错误的，因为先把0值去掉后产生的index和先排序再去0得到的index是不同的
%     NCC_matrix = (Mat_feat_2'*Mat_feat_1)/(window_length^2-1);
%     NCC_matrix(NCC_matrix < min_corr) = 0;
%     pre_sorted = reshape(NCC_matrix, no_pts2*no_pts1, 1); % size(Mat_feat_1) = (window_length^2, no_pts2)
%     pre_sorted(pre_sorted == 0) = [];
%     [sorted, sorted_index] = sort(pre_sorted(:), 1, 'descend'); %别忘了用:全部筛选，此时不要写dim很容易犯错 !!!

    NCC_matrix = (Mat_feat_2'*Mat_feat_1)/(window_length^2-1); % size(NCC_matrix) == (no_pts2, no_pts1);
    NCC_matrix(NCC_matrix < min_corr) = 0;

    [sorted, sorted_index] = sort(NCC_matrix(:), 'descend'); %别忘了用:全部筛选: !!!
    sorted_index(sorted == 0) = [];

    %% Correspondeces
    % 目标：通过sorted_index中的特征标号找到Ftp中存储的特征点坐标
    % sorted_index->NCC_matrix(no_pts2*no_pts1)->Mat_feat对应的坐标->Ftp对应的index


%     cor = zeros(4, size(sorted_index, 1));
%     count = 1;
%     for i = 1:size(sorted_index, 1)
%         NCC_index_product = sorted_index(i);
%         x_NCC_matrix = ceil(NCC_index_product/no_pts2);
%         y_NCC_matrix = ceil(NCC_index_product-no_pts2*(x_NCC_matrix-1));
% %         if all(Mat_feat_1(:, x_NCC_matrix) == 0) || all(Mat_feat_2(:, y_NCC_matrix) == 0)
% %             continue;
% %         end
%         if any(Mat_feat_1(:, x_NCC_matrix) ~= 0) && any(Mat_feat_2(:, y_NCC_matrix) ~= 0)
%             cor(:, count) = [Ftp1(1, x_NCC_matrix); Ftp1(2, x_NCC_matrix); Ftp2(1, y_NCC_matrix); Ftp2(2, y_NCC_matrix)];
%             count = count+1;
%             Mat_feat_1(:, x_NCC_matrix) = 0;
%             Mat_feat_2(:, y_NCC_matrix) = 0;
%         else
%             continue
%         end
%         % 当I1中的一对点通过最高的correspondence匹配到I2中的一对点时，就将该对点由norm_windows产生Mat_feat
%         % 中所对应的列向量删去，避免重复匹配 
% 
%     end
%     cor(:, all(cor==0, 1)) = [];

%     cor = zeros(4, size(sorted_index, 1));
%     for i = 1:size(sorted_index, 1)
%         NCC_index_product = sorted_index(i);
%         x_NCC_matrix = ceil(NCC_index_product/no_pts2);
%         y_NCC_matrix = ceil(NCC_index_product-no_pts2*(x_NCC_matrix-1));
%         if all(Mat_feat_1(:, x_NCC_matrix) == 0) || all(Mat_feat_2(:, y_NCC_matrix) == 0)
%             continue;
%         end
% 
%         cor(:, i) = [Ftp1(1, x_NCC_matrix); Ftp1(2, x_NCC_matrix); Ftp2(1, y_NCC_matrix); Ftp2(2, y_NCC_matrix)];
%         % 当I1中的一对点通过最高的correspondence匹配到I2中的一对点时，就将该对点由norm_windows产生Mat_feat
%         % 中所对应的列向量删去，避免重复匹配 
%         Mat_feat_1(:, x_NCC_matrix) = 0;
%         Mat_feat_2(:, y_NCC_matrix) = 0;
%     end
%     cor(:, all(cor==0, 1)) = [];


    cor = zeros(4, size(sorted_index, 1));
    for i = 1:size(sorted_index, 1)
        NCC_index_product = sorted_index(i);
        x_NCC_matrix = ceil(NCC_index_product/no_pts2);
        y_NCC_matrix = ceil(NCC_index_product-no_pts2*(x_NCC_matrix-1));
        if all(Mat_feat_1(:, x_NCC_matrix) == 0) || all(Mat_feat_2(:, y_NCC_matrix) == 0)
            continue;
        end

        cor(:, i) = [Ftp1(1, x_NCC_matrix); Ftp1(2, x_NCC_matrix); Ftp2(1, y_NCC_matrix); Ftp2(2, y_NCC_matrix)];
        % 当I1中的一对点通过最高的correspondence匹配到I2中的一对点时，就将该对点由norm_windows产生Mat_feat
        % 中所对应的列向量删去，避免重复匹配 
        Mat_feat_1(:, x_NCC_matrix) = 0;
        Mat_feat_2(:, y_NCC_matrix) = 0;
    end

    %% Visualize the correspoinding image point pairs
    if do_plot == true
        imshow(Im1/255);
        hold on
        imshow(Im2/255);
        hold on
        alpha(0.5);  % add transparancy of 50%
        for i=1:size(cor,2)
            plot(cor(1,:), cor(2,:), 'r*', 'MarkerSize', 6, 'LineWidth',1);
            plot(cor(3,:), cor(4,:), 'b*', 'MarkerSize', 6, 'LineWidth',1);
            plot([cor(1,i),cor(3,i)],[cor(2,i),cor(4,i)], 'g-', 'LineWidth',1);
            hold on
        end
    end
end