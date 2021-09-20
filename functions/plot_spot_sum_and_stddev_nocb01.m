function plot_spot_sum_and_stddev_nocb01(i, vals, dev_img, fnames, numfids)

%i = spot_nr;
name = fnames(i).name;

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
f1.Name = ['spot number ' num2str(i) ' of ' num2str(numfids) ' is ' name ];

subplot(1, 2, 1)
       
    imshow(vals{i},[min(min(vals{i})) max(max(vals{i}))], 'InitialMagnification', 'fit'); hold on
    title('sum image');
    
subplot(1,2,2);
    imshow(dev_img(:,:,i),[min(min(min(dev_img(:,:,i)))) max(max(max(dev_img(:,:,i))))], 'InitialMagnification', 'fit'); hold on
    title('standard deviation');
    
hold off