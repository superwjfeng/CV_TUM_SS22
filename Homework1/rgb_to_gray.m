function gray_image = rgb_to_gray(input_image)
    % This function is supposed to convert a RGB-image to a grayscale image.
    % If the image is already a grayscale image directly return it.

    % 以这种错误写是错误的，但输入的是灰度图像时，input_iamge只有2个Channels
    % 读不到第三个Channels的信息
%     [r, c, channels] = size(input_image)
%     if channels == 1
%         gray_image = input_iamge;
%         return;
%     end
%     R = double(input_image(:,:,1));
%     G = double(input_image(:,:,2));
%     B = double(input_image(:,:,3));
%     gray_image = uint8(0.299*R + 0.587*G + 0.114*B);

    [ ~, ~, channels] = size(input_image);
    if channels == 3
        R = double(input_image(:,:,1));
        G = double(input_image(:,:,2));
        B = double(input_image(:,:,3));
        gray_image = uint8(0.299*R + 0.587*G + 0.114*B);
    else
        gray_image = input_image;
    end
end
