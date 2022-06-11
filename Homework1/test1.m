img = imread('test.png');
gray = rgb_to_gray(img);
features = harris_detector(gray, 'segment_length', 9, 'k', 0.06, 'do_plot', true);