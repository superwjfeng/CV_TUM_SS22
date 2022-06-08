function [min_dist, tile_size, N] = harris_detector_upgrade(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.

    %% Input parser
    p = inputParser;
    addRequired(p, 'input_image');
    
    default_segment_length = 15;
    addOptional(p, 'segment_length', default_segment_length, @(x) isnumeric(x) && (x > 1) && (mod(x,2)~=0));

    default_k = 0.05;
    addOptional(p, 'k', default_k, @(x) isnumeric(x) && x >= 0 && x <= 1);

    default_tau = 10^6;
    addOptional(p, 'tau', default_tau, @(x) isnumeric(x) && x > 0);

    default_do_plot = false;
    addOptional(p, 'do_plot',default_do_plot, @(x) islogical(x));

    default_min_dist = 20;
    addOptional(p, 'min_dist', default_min_dist, @(x) isnumeric(x) && x>=1);

    default_tile_size = [200, 200];
    addOptional(p, 'tile_size', default_tile_size, @(x) isnumeric(x));

    default_N = 5;
    addOptional(p, 'N', default_N, @(x) isnumeric(x) && x>=1);

    parse(p,'input_image', varargin{:});

    if length(p.Results.tile_size) == 2
        Results.tile_size = p.Results.tile_size;
    else
        height = p.Results.tile_size;
        Results.tile_size = [height, 200];
    end

    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist = p.Results.min_dist;
    tile_size = Results.tile_size;
    N = p.Results.N;

    %%
    % 总体思路：角点是通过x和y方向都存在较大变化（即x和y方向都是边缘）找出来的
    % 1、通过sobel卷积核对图像作平滑+一阶导（边缘提取）得到Ix和Iy，即找出可能的边
    % 2、计算二阶矩矩阵，通过高斯窗函数和二阶矩的卷积得到。光有边是不够的，要得到角点必须要知道每个像素点在临近的窗宽（选定）
    % 中一阶导的变化关系，即求x的移动值与其一阶导变化量的关系的损失函数E，如果损失函数为0，则说明至少有一边是没有变化的
    % 每个像素的E损失矩阵在泰勒分解后可以将其写成与黑塞矩阵的关联式
    % 3、计算Harris判别式，当黑塞正定时（正定等价于该矩阵的两个特征值大于0）可以知道E不等于0，但由于分别判断lambda1
    % lambda2比较麻烦，因此给出经验公式R，即Harris判别式
    % 4、Threshold R，很多R比较小时正定性不佳
    % 5、对R做非最大化抑制，选出最优
    % 最终将问题转换为了求判别式
    %% Preparation for feature extraction
    % Check if it is a grayscale image
    [r, c, ch] = size(input_image);
    if ch ~= 1
        error('Image format has to be NxMx1');
    end

    % Approximation of the image gradient
    im2double(input_image);
    [Ix, Iy] = sobel_xy(input_image);
    
    % Weighting
    % 设置Gauss窗函数，调整每个点对总的损失的影响，突出中心
    w = fspecial('gaussian',[1 segment_length], segment_length/5);

    % Harris Matrix G(The summation Hessian matrix and windows function of every pixel)
    Gxx = conv2(w, w, Ix.^2, 'same'); % 拆成两个小卷积核
    Gyy = conv2(w, w, Iy.^2, 'same');
    Gxy = conv2(w, w, Ix.*Iy, 'same');

    %% Feature extraction with the Harris measurement
    H = (Gxx.*Gyy - Gxy.^2) - k*(Gxx+Gyy).^2;
    mask = zeros(size(H)); 
    mask((ceil(segment_length/2)+1):(size(H,1)-ceil(segment_length/2)), (ceil(segment_length/2)):(size(H,2)-ceil(segment_length/2))) = 1;
    corners = H.*mask; %去除边缘进行零填充后卷积的性质不好的点
    
    corners(corners<tau) = 0;
    [row, col] = find(corners);
    % features = [col'; row']; % col为x，row为y 
    features = [col, row]';
     
    %% Plot
%     if do_plot == true
%         figure('Name', 'Harris Plot');
%         imshow(input_image);
%         hold on;
%         for i = 1:size(features, 2) % 被误导了，plot和scatter不需要使用循环，直接输入x和y向量后画一次就出来了
%             scatter(features(1, i), features(2, i), 'o');  
%         end
%     end

    if do_plot == true
        imshow(input_image);
        hold on;
        plot(features(1,:), features(2,:), 'ro');
    end

end