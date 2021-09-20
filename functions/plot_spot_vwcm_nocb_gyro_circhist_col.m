function [class, dwell, pick, Theta_D_shifted, rotation] = plot_spot_vwcm_nocb_gyro_circhist_col(spot_nr, K, data, fnamesfits, dev_img)
%PLOT_SPOT Plots one spot with automatic clustering in three sets

rad1 = 10;
rad2 = 10;
i = spot_nr;
xydata = data{i}{1,1}.vwcm.pos;

f1 = figure('units','normalized','outerposition',[0 0 1 1]);


    class = kmeans([xydata(:,1),xydata(:,2)], 3);
    plot(xydata(class==1,1),xydata(class==1,2),'b.'), hold on
    scatter_kde(xydata(:,1),(xydata(:,2)), 'filled', 'MarkerSize', 8);
    %hold on
    
    %Plots mean dot of class 1
    class1_mean = [mean(xydata(class==1,1)), mean(xydata(class==1,2))];
    plot(class1_mean(1), class1_mean(2), 'xk')

    %plot(xydata(class==2,1),xydata(class==2,2),'r.')
    class2_mean = [mean(xydata(class==2,1)), mean(xydata(class==2,2))];
   
    %Plots mean dot of class 2
    plot(class2_mean(1), class2_mean(2), 'xk')
    
    %plot(xydata(class==3,1),xydata(class==3,2),'g.')
    class3_mean = [mean(xydata(class==3,1)), mean(xydata(class==3,2))];
    
    %Plots mean dot of class 3
    plot(class3_mean(1), class3_mean(2), 'xk')

    % Plots channel mean
    ch1_mean = mean(xydata);
    plot(ch1_mean(1), ch1_mean(2), 'ok','MarkerSize',15)
        
    % Set axis
    %mean_current = mean(data{mov_nr}{i,1}.vwcm.pos);
    mean_current = [xydata(2,1), xydata(2,2)];
    xlim ([mean_current(1)-rad1,mean_current(1)+rad1])
    ylim ([mean_current(2)-rad1,mean_current(2)+rad1])
    axis square
    
    %Write title
    %d(i) = sqrt((class1_mean(1)-class2_mean(1)).^2+(class1_mean(2)-class2_mean(2)).^2);
    title(['Red channel spot ' num2str(i) ' of ' num2str(K) ' in ' num2str(fnamesfits(i).name(1:end-8))]  );
    
    %  Takes three inputs for class1 mean, class2 mean, class3 mean
    [x,y] = ginput(3);
    pick = [x,y];
    
    new_mean = mean(pick);

    close gcf
    
    
f2 = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3, 8, [1:2 9:10])
       
    %class = kmeans([data{mov_nr}{i,1}.vwcm.pos(:,1),data{mov_nr}{i,1}.vwcm.pos(:,2)], 2);
    class = kmeans([xydata(:,1),xydata(:,2)], 3, 'MaxIter', 1, 'Start', pick);
    plot(xydata(class==1,1),xydata(class==1,2),'b.'), hold on 
    plot(xydata(class==2,1),xydata(class==2,2),'r.')
    plot(xydata(class==3,1),xydata(class==3,2),'g.')
    class3_mean = [mean(xydata(class==3,1)), mean(xydata(class==3,2))];

   


    % Plots channel mean
    new_mean = mean(pick);
    plot(new_mean(1), new_mean(2), 'bo','MarkerSize',5)
    plot([pick(1,1) pick(2,1) pick(3,1)], [pick(1,2) pick(2,2) pick(3,2)],'bx','LineWidth',2)
    
    % Set axis
    %mean_current = mean(data{mov_nr}{i,1}.vwcm.pos);
    %mean_current = [ data{mov_nr}{i,1}.vwcm.pos(1,1), data{mov_nr}{i,1}.vwcm.pos(1,2)];
    xlim ([new_mean(1)-rad2,new_mean(1)+rad2])
    ylim ([new_mean(2)-rad2,new_mean(2)+rad2])
    set (gca, 'Ydir', 'reverse')
    axis square
    
    %Write title
    % d(i) = sqrt((class1_mean(1)-class2_mean(1)).^2+(class1_mean(2)-class2_mean(2)).^2);
    %title(['Red channel spot ' num2str(i) ' of ' num2str(size(data{mov_nr},1)) ' in movie ' num2str(mov_nr) ' distance = ' num2str(d(i)) ' px']  );
    title(['Red channel spot ' num2str(i) ' of ' num2str(K) ' in movie ' num2str(fnamesfits(i).name(1:end-8))]  );

    
subplot(3, 8, [3:4 11:12])
   
    imshow(dev_img(:,:,i),[min(min(min(dev_img(:,:,i)))) max(max(max(dev_img(:,:,i))))], 'InitialMagnification', 'fit'); hold on
    title('standard deviation');

    
%subplot(3, 8, [5:8])
subplot (3, 8, [13:16])
for j = [1:length(data{i}{1,1}.vwcm.pos(:,1))]
 
        if (data{i}{1,1}.vwcm.pos(j,1) == 0)
            d(j) = 0;
        else
            d(j) = sqrt((data{i}{1,1}.vwcm.pos(j,1)-new_mean(1)).^2+(data{i}{1,1}.vwcm.pos(j,2)-new_mean(2)).^2);
        end
end

plot(d)
mean_d = mean(d(d~=0));
hline(2.76)
ylim ([0,10])
title(['Normalized distance from mean value, mean is ' num2str(mean_d)]);


%calculate dwell times for 3 different classes
f1 = find(diff([false, class'==1, false])~=0);
L1 = f1(2:2:end)-f1(1:2:end-1);

f2 = find(diff([false, class'==2, false])~=0);
L2 = f2(2:2:end)-f2(1:2:end-1);

f3 = find(diff([false, class'==3, false])~=0);
L3 = f3(2:2:end)-f3(1:2:end-1);


subplot(3, 8, 5);
    scatter_kde(xydata(:,1),(xydata(:,2)*-1), 'filled', 'MarkerSize', 8);
    xlim ([new_mean(1)-5,new_mean(1)+5])
    ylim ([new_mean(2)*-1-5,new_mean(2)*-1+5])
    axis square

subplot(3, 8, 6);
hist(L1, 50);
title('cluster 1, blue');
xlabel('dwell time (#frames)')

subplot(3, 8, 7);
hist(L2, 50);
title('cluster 2, red');
xlabel('dwell time (#frames)')

subplot(3, 8, 8);
hist(L3, 50);
title('cluster 3, green');
xlabel('dwell time (#frames)')


% subplot(3, 8, [13:16])
% plot(data{i}{1,1}.med_itrace)
% title('Intensity trace');

%dwell = {[L1, zeros(1, 1000-length(L1))], [L2, zeros(1, 1000-length(L2))], [L3, zeros(1, 1000-length(L3))]};
dwell = {[L1, zeros(1, 10000-length(L1))], [L2, zeros(1, 10000-length(L2))], [L3, zeros(1, 10000-length(L3))]};


rotation = diff(class);
rotation(rotation==-2)=1;
rotation(rotation==2)=-1;

Theta_cum(1) = 0;

 normal = pick(1,:) - new_mean;
 normal = normal./norm(normal);
 [Theta_norm,Rho] = cart2pol(normal(1,1),normal(1,2));
 Theta_norm = wrapTo2Pi(Theta_norm);
 Theta_norm_shifted = Theta_norm + (60*pi/180);
 %Theta_norm = 0;
 Theta_norm_D = Theta_norm*180/pi;
 Theta_norm_D_shifted = Theta_norm_shifted*180/pi;
 R = [cosd(Theta_norm_D) -sind(Theta_norm_D); sind(Theta_norm_D) cosd(Theta_norm_D)];
 R_shifted = [cosd(Theta_norm_D_shifted) -sind(Theta_norm_D_shifted); sind(Theta_norm_D_shifted) cosd(Theta_norm_D_shifted)];


for k = 1:(size(data{i}{1,1}.vwcm.pos,1)-1) % Goes through all frames
    
  vector = data{i}{1,1}.vwcm.pos(k,:) - new_mean;
  vector = vector./norm(vector);
  vector = vector*R;
  vector_shifted = vector*R_shifted;
  
  vector2 = data{i}{1,1}.vwcm.pos(k+1,:) - new_mean;
  vector2 = vector2./norm(vector2);
  vector2 = vector2*R;
  vector2_shifted = vector2*R_shifted;
  
  
  %CosTheta(k) = dot(normal,vector)/(norm(normal)*norm(vector));
  %Theta_D(k) = atan(vector(1,2)./vector(1,1))*180/pi;
  
  [Theta(k),Rho] = cart2pol(vector(1,1),vector(1,2));
  [Theta_shifted(k),Rho] = cart2pol(vector_shifted(1,1),vector_shifted(1,2));

  %Theta(k) = Theta(k) + pi;
  
  %Theta_bet(k) = dot(vector2,vector)/(norm(vector2)*norm(vector));
  %Theta_2_D(k) = atan(vector2(1,2)./vector2(1,1))*180/pi;
  
  [Theta_2(k),Rho] = cart2pol(vector2(1,1),vector2(1,2));
  %Theta_2(k) = Theta_2(k) + pi;
  
  Theta(k) = wrapTo2Pi(Theta(k));
  Theta_shifted(k) = wrapTo2Pi(Theta_shifted(k));
  Theta_2(k) = wrapTo2Pi(Theta_2(k));
  
  %Theta inbetween
  Theta_bet(k) = Theta_2(k) - Theta(k);
  
  %ensure shortest distence has been chosen (i.e. forbid theta_bet>180�)
  if abs(Theta_bet(k)) < pi
     Theta_bet(k) = Theta_bet(k);
  else
     if Theta_bet(k) < 0
         Theta_bet(k) = Theta_bet(k) + 2*pi;
     else
         Theta_bet(k) = Theta_bet(k) - 2*pi;
     end
  end
 
  Theta_cum(k+1) = Theta_bet(k) + Theta_cum(k);
  
end

  Theta_D = Theta*180/pi;
  Theta_D_shifted = Theta_shifted*180/pi;
  Theta_2_D = Theta_2*180/pi;
  Theta_bet_D = Theta_bet*180/pi;
  Theta_cum_D = Theta_cum*180/pi;
  
% Plot Theta and Histogram
subplot(3, 8, [17:20])

    %medcostheta1 = medfilt1(Theta_D,5);
    hist(Theta_D_shifted, 360);
    ylim ([0,80])
    title('radial distribution of Theta');

subplot(3, 8, [21 24])

    %medcostheta2 = medfilt1(Theta_cum,5);
    framenr = [1:length(Theta_cum_D)];
    plot(Theta_cum_D, 'black'), hold on;
    plot(framenr(class==1), Theta_cum_D(class==1),'b.');
    plot(framenr(class==2), Theta_cum_D(class==2),'r.');
    plot(framenr(class==3), Theta_cum_D(class==3),'g.');
    title('Theta cumulative');

hold off