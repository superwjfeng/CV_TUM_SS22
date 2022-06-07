function [Fx, Fy] = sobel_xy(input_image)
    % In this function you have to implement a Sobel filter 
    % that calculates the image gradient in x- and y- direction of a grayscale image.
    sobel_xkernel = [-1 0 1;-2 0 2;-1 0 1];
    sobel_xkernel_flip = flip(sobel_xkernel, 2); %图像处理中卷积核要进行翻转
    sobel_ykernel = [1 2 1;0 0 0;-1 -2 -1];
    sobel_ykernel_flip = flip(sobel_ykernel, 2);
    Fx = conv2(input_image, sobel_xkernel_flip, 'same');
    % flip(A)默认第二位参数为1，此时按行flip，当参数为2按行flip
    Fy = conv2(input_image, sobel_ykernel_flip, 'same');
end