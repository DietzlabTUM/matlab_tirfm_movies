function [class, dwell, pick, Theta_D_shifted, Theta_cum_D, Theta_bet_D, rotation, net_rot, ind, turns, new_mean] = plot_spot_eField_nocb_gyro_col(spot_nr, K, data, eField, varargin)
%PLOT_SPOT Plots one spot with automatic clustering in three sets

rad1 = 10;
rad2 = 10;
i = spot_nr;
xydata = data{i}.fit;
xydata_nonan = xydata((all((~isnan(xydata)),2)),:);
pick = [0 0];
movie_length = length(xydata);

%%   %calculate angle of lever arm
    
    if length(pick(:,1)) == 1
        [new_mean(1) new_mean(2) R] = circle_fit(xydata_nonan(:,1), xydata_nonan(:,2));
        %new_mean = nanmean(xydata);
    else
        new_mean = mean(pick);
    end
    
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
    
    for k = 1:(size(xydata,1)-1) % Goes through all frames
        
        vector = xydata(k,:) - new_mean;
        vector = vector./norm(vector);
        vector = vector*R;
        vector_shifted = vector*R_shifted;
        
        vector2 = xydata(k+1,:) - new_mean;
        vector2 = vector2./norm(vector2);
        vector2 = vector2*R;
        vector2_shifted = vector2*R_shifted;
        
        [Theta(k),Rho] = cart2pol(vector(1,1),vector(1,2));
        [Theta_shifted(k),Rho] = cart2pol(vector_shifted(1,1),vector_shifted(1,2));
        
        [Theta_2(k),Rho] = cart2pol(vector2(1,1),vector2(1,2));
        
        Theta(k) = wrapTo2Pi(Theta(k));
        Theta_shifted(k) = wrapTo2Pi(Theta_shifted(k));
        Theta_2(k) = wrapTo2Pi(Theta_2(k));
        
        %Theta inbetween
        Theta_bet(k) = Theta_2(k) - Theta(k);
        
        %ensure shortest distence has been chosen (i.e. forbid theta_bet>180ï¿½)
        if abs(Theta_bet(k)) < pi
            Theta_bet(k) = Theta_bet(k);
        else
            if Theta_bet(k) < 0
                Theta_bet(k) = Theta_bet(k) + 2*pi;
            else
                Theta_bet(k) = Theta_bet(k) - 2*pi;
            end
        end
        
        if isnan(Theta_bet(k))
            Theta_cum(k+1) = Theta_cum(k);
        else
            Theta_cum(k+1) = Theta_bet(k) + Theta_cum(k);
        end
        
    end
    
    Theta_D = Theta*180/pi;
    Theta_D_shifted = Theta_shifted*180/pi;
    Theta_2_D = Theta_2*180/pi;
    Theta_bet_D = Theta_bet*180/pi;
    Theta_cum_D = Theta_cum*180/pi;
    
    %get angle of picks
    for k = 1:length(pick(:,1)) % Goes through all picks
        
        vector_p = pick(k,:) - new_mean;
        vector_p = vector_p./norm(vector_p);
        vector_p = vector_p*R;
        vector_p_shifted = vector_p*R_shifted;
        
        [Theta_pick(k),Rho] = cart2pol(vector_p(1,1),vector_p(1,2));
        [Theta_pick_shifted(k),Rho] = cart2pol(vector_p_shifted(1,1),vector_p_shifted(1,2));
        
        Theta_pick(k) = wrapTo2Pi(Theta_pick(k));
        Theta_pick_shifted(k) = wrapTo2Pi(Theta_pick_shifted(k));
        
        
    end
    
    Theta_pick_D = Theta_pick*180/pi;
    Theta_pick_D_shifted = Theta_pick_shifted*180/pi;
    
  
%% find max
p = length(pick(:,1));
    if p == 6
        %find 6 maxima in theta_hist and assign clusters
        h = hist(Theta_D_shifted, 360);
        s = smooth(h);
        s_spec(1:45) = 0;
        s_spec(46:405) = s;
        s_spec(406:450) = 0;
        m_spec = islocalmax(s_spec,'MaxNumExtrema',6, 'MinSeparation', 45);
        m = m_spec(46:405);
        w = [0:360];
        maxima = w(m);
        max_values = h(m);
        
        if 360-maxima(end)+maxima(1) < 45
            if max_values(1) > max_values(end)
                maxima = maxima(1:end-1);
            else
                maxima = maxima(2:end);
            end
        end
        
        if maxima(1) == 0
            maxima(1) = 1;
        end
        m_new = zeros(1,360);
        m_new(1,maxima) = 1;
        m_new = logical(m_new);

        start = zeros(length(maxima),2);
        ind=0;
        for n = 1:length(maxima)
            [M,I] = min(abs(Theta_D_shifted-maxima(n)));
            ind(n) = I;
            start(n, :) = xydata(I, :);
        end

        number_of_max = length(ind);

        % replace missing max with picks
        if number_of_max ~= p
        miss = ismembertol(Theta_pick_D_shifted, maxima, 0.025);
        miss_pos = find(miss == 0);
            if length(miss_pos) + number_of_max == p
                start = [start; pick(miss_pos, :)];
            else
                last_max_repl = maxima(end)-360;
                max_repl = maxima;
                max_repl(end) = abs(last_max_repl);
                miss_repl = ismembertol(Theta_pick_D_shifted, max_repl, 0.025);
                miss_repl_pos = find(miss_repl == 0);
                start = [start; pick(miss_repl_pos, :)];
                if length(start) ~= p
                    start = pick;
                end
            end
        end 
        
    elseif p == 1
        h = hist(Theta_D_shifted, 360);
        s = smooth(h);
        s_spec(1:45) = 0;
        s_spec(46:405) = s;
        s_spec(406:450) = 0;
        m_spec = islocalmax(s_spec,'MaxNumExtrema', 6, 'MinSeparation', 45);
        m = m_spec(46:405);
        w = [0:360];
        maxima = w(m);
        max_values = h(m);
        
        
        if 360-maxima(end)+maxima(1) < 45
            if max_values(1) > max_values(end)
                maxima = maxima(1:end-1);
            else
                maxima = maxima(2:end);
            end
        end
        
        if maxima(1) == 0
            maxima(1) = 1;
        end
        m_new = zeros(1,360);
        m_new(1,maxima) = 1;
        m_new = logical(m_new);
        
        start = zeros(length(maxima),2);
        ind=0;
        for n = 1:length(maxima)
            [M,I] = min(abs(Theta_D_shifted-maxima(n)));
            ind(n) = I;
            start(n, :) = xydata(I, :);
        end

        number_of_max = length(ind)
        
    else
        
        %find p maxima in theta_hist and assign clusters
        h = hist(Theta_D_shifted, 360);
        s = smooth(h);
        s_spec(1:45) = 0;
        s_spec(46:405) = s;
        s_spec(406:450) = 0;
        m_spec = islocalmax(s_spec,'MaxNumExtrema',p, 'MinSeparation', 360/p*0.75);
        m = m_spec(46:405);
        w = [0:360];
        maxima = w(m);
        max_values = h(m);
        
        if 360-maxima(end)+maxima(1) < 360/p*0.75
            if max_values(1) > max_values(end)
                maxima = maxima(1:end-1);
            else
                maxima = maxima(2:end);
            end
        end
        
        if maxima(1) == 0
            maxima(1) = 1;
        end
        m_new = zeros(1,360);
        m_new(1,maxima) = 1;
        m_new = logical(m_new);
        
        start = zeros(length(maxima),2);
        ind=0;
        for n = 1:length(maxima)
            [M,I] = min(abs(Theta_D_shifted-maxima(n)));
            ind(n) = I;
            start(n, :) = xydata(I, :);
        end

        number_of_max = length(ind);
        
        if number_of_max ~= p
        miss = ismembertol(Theta_pick_D_shifted, maxima, 0.025);
        miss_pos = find(miss == 0);
            if length(miss_pos) + number_of_max == p
                start = [start; pick(miss_pos, :)];
            else
                last_max_repl = maxima(end)-360;
                max_repl = maxima;
                max_repl(end) = abs(last_max_repl);
                miss_repl = ismembertol(Theta_pick_D_shifted, max_repl, 0.025);
                miss_repl_pos = find(miss_repl == 0);
                start = [start; pick(miss_repl_pos, :)];
                if length(start) ~= p
                    start = pick;
                end
            end
        end     
    end
    
%%    %fit 6 clusters according to manual selection of cluster centers
    if p ~= 1
        class = kmeans([xydata(:,1),xydata(:,2)], p, 'Start', start);
    else
        
        class = kmeans([xydata(:,1),xydata(:,2)], number_of_max, 'Start', start);
    end
    
    %sort classes with highest occupancy
    class_tab = tabulate(class);
    class_sort = sortrows(class_tab,3, 'descend');
    
    
    %calculate dwell times for 6 different classes and add to dwell
    f1 = find(diff([false, class'==1, false])~=0);
    L1 = f1(2:2:end)-f1(1:2:end-1);

    f2 = find(diff([false, class'==2, false])~=0);
    L2 = f2(2:2:end)-f2(1:2:end-1);

    f3 = find(diff([false, class'==3, false])~=0);
    L3 = f3(2:2:end)-f3(1:2:end-1);
    
    f4 = find(diff([false, class'==4, false])~=0);
    L4 = f4(2:2:end)-f4(1:2:end-1);
    
    f5 = find(diff([false, class'==5, false])~=0);
    L5 = f5(2:2:end)-f5(1:2:end-1);
    
    f6 = find(diff([false, class'==6, false])~=0);
    L6 = f6(2:2:end)-f6(1:2:end-1);

    dwell = {[L1, zeros(1, 8000-length(L1))], [L2, zeros(1, 8000-length(L2))], [L3, zeros(1, 8000-length(L3))], [L4, zeros(1, 8000-length(L4))], [L5, zeros(1, 8000-length(L5))], [L6, zeros(1, 8000-length(L6))]};

    %determine counter-clockwise (-1) and clockwise (1) switches
    rotation = diff(class);
    rotation(rotation==-2)=1;
    rotation(rotation==2)=-1;
    rotation_corr=rotation*(-1);
    net_rot = nansum(rotation_corr);
    turns = Theta_cum_D(end)/360;
    
    %percentage of fitted frames
    nonan = nnz(~isnan(data{i}.fit(:,1)));
    perc = nonan/length(data{i}.fit)*100;
    
%%    % eField plot
    
    % static with wait time
    % wait = eField(1); dur = eField(2); iter = eField(3); tpf = eField(4),
    % vol = eField(5)
    if length(eField) == 5
        waitframes = 100*eField(1)*1000/eField(4);
        period = eField(5)*[ones(1, 100*eField(2)/eField(4)), -ones(1, 100*eField(2)/eField(4))];
        fullperiod = repmat(period, [1, eField(3)]);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:waitframes) = 0;
        ef(1, waitframes+1:waitframes+length(fullperiod)) = fullperiod;
        ef(1, waitframes+length(fullperiod)+1:end)=0;
    
    % ramp without wait inbetween
    % wait = eField(1); dur1 = eField(2); dur2 = eField(3); dur3 = eField(4);
    % dur4 = eField(5); numframes = eField(6); tpf = eField(7); vol = eField(8)
    elseif length(eField) == 8
        waitframes = 100*eField(1)*1000/eField(7);
        period1 = eField(8)*[ones(1, 100*eField(2)/eField(7)), -ones(1, 100*eField(2)/eField(7))];
        period2 = eField(8)*[ones(1, 100*eField(3)/eField(7)), -ones(1, 100*eField(3)/eField(7))];
        period3 = eField(8)*[ones(1, 100*eField(4)/eField(7)), -ones(1, 100*eField(4)/eField(7))];
        period4 = eField(8)*[ones(1, 100*eField(5)/eField(7)), -ones(1, 100*eField(5)/eField(7))];
        fullperiod1 = repmat(period1, [1, (eField(6)/(eField(2)/eField(7)*2))]);
        fullperiod2 = repmat(period2, [1, (eField(6)/(eField(3)/eField(7)*2))]);
        fullperiod3 = repmat(period3, [1, (eField(6)/(eField(4)/eField(7)*2))]);
        fullperiod4 = repmat(period4, [1, (eField(6)/(eField(5)/eField(7)*2))-1]);
            ef = zeros(1, 100*length(Theta_cum_D));
            ef(1, 1:waitframes) = 0;
            ef(1, waitframes+1:waitframes+100*eField(6)) = fullperiod1;
            ef(1, waitframes+100*eField(6)+1:waitframes+200*eField(6)) = fullperiod2;
            ef(1, waitframes+200*eField(6)+1:waitframes+300*eField(6)) = fullperiod3;
            ef(1, waitframes+300*eField(6)+1:end) = fullperiod4;
            
    % long ramp without wait inbetween (freq)
    % wait = eField(1); dur1 = eField(2); dur2 = eField(3); dur3 = eField(4);
    % dur4 = eField(5); dur5 = eField(6); dur6 = eField(7); numframes = eField(8);
    % tpf = eField(9); vol = eField(10)
    elseif length(eField) == 12
        waitframes = 100*eField(1)*1000/eField(9);
        period1 = eField(10)*[ones(1, 100*eField(2)/eField(9)), -ones(1, 100*eField(2)/eField(9))];
        period2 = eField(10)*[ones(1, 100*eField(3)/eField(9)), -ones(1, 100*eField(3)/eField(9))];
        period3 = eField(10)*[ones(1, 100*eField(4)/eField(9)), -ones(1, 100*eField(4)/eField(9))];
        period4 = eField(10)*[ones(1, 100*eField(5)/eField(9)), -ones(1, 100*eField(5)/eField(9))];
        period5 = eField(10)*[ones(1, 100*eField(6)/eField(9)), -ones(1, 100*eField(6)/eField(9))];
        period6 = eField(10)*[ones(1, 100*eField(7)/eField(9)), -ones(1, 100*eField(7)/eField(9))];
        fullperiod1 = repmat(period1, [1, (eField(8)/(eField(2)/eField(9)*2))]);
        fullperiod2 = repmat(period2, [1, (eField(8)/(eField(3)/eField(9)*2))]);
        fullperiod3 = repmat(period3, [1, (eField(8)/(eField(4)/eField(9)*2))]);
        fullperiod4 = repmat(period4, [1, (eField(8)/(eField(5)/eField(9)*2))]);
        fullperiod5 = repmat(period5, [1, (eField(8)/(eField(6)/eField(9)*2))]);
        fullperiod6 = repmat(period6, [1, (eField(8)/(eField(7)/eField(9)*2))]);
            ef = zeros(1, 100*length(Theta_cum_D));
            ef(1, 1:waitframes) = 0;
            ef(1, waitframes+1:waitframes+100*eField(8)) = fullperiod1;
            ef(1, waitframes+100*eField(8)+1:waitframes+200*eField(8)) = fullperiod2;
            ef(1, waitframes+200*eField(8)+1:waitframes+300*eField(8)) = fullperiod3;
            ef(1, waitframes+300*eField(8)+1:waitframes+400*eField(8)) = fullperiod4;
            ef(1, waitframes+400*eField(8)+1:waitframes+500*eField(8)) = fullperiod5;
            ef(1, waitframes+500*eField(8)+1:waitframes+600*eField(8)) = fullperiod6;
            
    % on-off with wait inbetween
    % waitbefore = eField(1); dur = eField(2); iter = eField(3); waitbetween = eField (4)
    % tpf = eField(5), vol = eField(6), mod = eField(7)
    elseif length(eField) == 7
        waitframes = 100*eField(1)*1000/eField(5);
        period = eField(6)*[ones(1, 100*eField(2)/eField(5)), -ones(1, 100*eField(2)/eField(5))];
        fullperiod = repmat(period, [1, eField(3)]);
        waitbetween = 100*eField(4)*1000/eField(5);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:waitframes) = 0;
        for q = 1:eField(7)
        ef(1, waitframes+(q-1)*length(fullperiod)+(q-1)*waitbetween+1:waitframes+q*length(fullperiod)+(q-1)*waitbetween) = fullperiod;
        ef(1, waitframes+q*length(fullperiod)+(q-1)*waitbetween+1:waitframes+q*length(fullperiod)+q*waitbetween) = 0;
        end
    
    % ramp with wait inbetween
    % waitbefore = eField(1); dur1 = eField(2); dur2 = eField(3); dur3 = eField(4);
    % vol1 = eField(5); vol2 = eField(6); vol3 = eField(7); numofframes = eField(8);
    % waitbetween = eField (9); tpf = eField(10)
    elseif length(eField) == 10
        waitframes = 100*eField(1)*1000/eField(10);
        period1 = eField(5)*[ones(1, 100*eField(2)/eField(10)), -ones(1, 100*eField(2)/eField(10))];
        period2 = eField(6)*[ones(1, 100*eField(3)/eField(10)), -ones(1, 100*eField(3)/eField(10))];
        period3 = eField(7)*[ones(1, 100*eField(4)/eField(10)), -ones(1, 100*eField(4)/eField(10))];
        fullperiod1 = repmat(period1, [1, (eField(8)/(eField(2)/eField(10)*2))]);
        fullperiod2 = repmat(period2, [1, (eField(8)/(eField(3)/eField(10)*2))]);
        fullperiod3 = repmat(period3, [1, (eField(8)/(eField(4)/eField(10)*2))]);
        waitbetween = 100*eField(9)*1000/eField(10);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:waitframes) = 0;
        ef(1, waitframes+1:waitframes+100*eField(8)) = fullperiod1;
        ef(1, waitframes+100*eField(8)+1:waitframes+100*eField(8)+waitbetween) = 0;
        ef(1, waitframes+100*eField(8)+waitbetween+1:waitframes+200*eField(8)+waitbetween) = fullperiod2;
        ef(1, waitframes+200*eField(8)+waitbetween+1:waitframes+200*eField(8)+2*waitbetween) = 0;
        ef(1, waitframes+200*eField(8)+2*waitbetween+1:waitframes+300*eField(8)+2*waitbetween) = fullperiod3;
        ef(1, waitframes+300*eField(8)+2*waitbetween+1:end) = 0;
        
    % ramp long with wait inbetween
    % waitbefore = eField(1); dur = eField(2); vol1 = eField(3); vol2 = eField(4);
    % vol3 = eField(5); vol4 = eField(6); vol5 = eField(7); vol6 = eField(8);
    % iter = eField(9); waitbetween = eField (10); tpf = eField(11)
    elseif length(eField) == 11
        waitframes = 100*eField(1)*1000/eField(11);
        period1 = eField(3)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        period2 = eField(4)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        period3 = eField(5)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        period4 = eField(6)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        period5 = eField(7)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        period6 = eField(8)*[ones(1, 100*eField(2)/eField(11)), -ones(1, 100*eField(2)/eField(11))];
        fullperiod1 = repmat(period1, [1, eField(9)]);
        fullperiod2 = repmat(period2, [1, eField(9)]);
        fullperiod3 = repmat(period3, [1, eField(9)]);
        fullperiod4 = repmat(period4, [1, eField(9)]);
        fullperiod5 = repmat(period5, [1, eField(9)]);
        fullperiod6 = repmat(period6, [1, eField(9)]);
        waitbetween = 100*eField(10)*1000/eField(11);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:waitframes) = 0;
        ef(1, waitframes+1:waitframes+length(fullperiod1)) = fullperiod1;
        ef(1, waitframes+length(fullperiod1)+1:waitframes+length(fullperiod1)+waitbetween) = 0;
        ef(1, waitframes+length(fullperiod1)+waitbetween+1:waitframes+2*length(fullperiod1)+waitbetween) = fullperiod2;
        ef(1, waitframes+2*length(fullperiod1)+waitbetween+1:waitframes+2*length(fullperiod1)+2*waitbetween) = 0;
        ef(1, waitframes+2*length(fullperiod1)+2*waitbetween+1:waitframes+3*length(fullperiod1)+2*waitbetween) = fullperiod3;
        ef(1, waitframes+3*length(fullperiod1)+2*waitbetween+1:waitframes+3*length(fullperiod1)+3*waitbetween) = 0;
        ef(1, waitframes+3*length(fullperiod1)+3*waitbetween+1:waitframes+4*length(fullperiod1)+3*waitbetween) = fullperiod4;
        ef(1, waitframes+4*length(fullperiod1)+3*waitbetween+1:waitframes+4*length(fullperiod1)+4*waitbetween) = 0;
        ef(1, waitframes+4*length(fullperiod1)+4*waitbetween+1:waitframes+5*length(fullperiod1)+4*waitbetween) = fullperiod5;
        ef(1, waitframes+5*length(fullperiod1)+4*waitbetween+1:waitframes+5*length(fullperiod1)+5*waitbetween) = 0;
        ef(1, waitframes+5*length(fullperiod1)+5*waitbetween+1:waitframes+6*length(fullperiod1)+5*waitbetween) = fullperiod6;
        ef(1, waitframes+6*length(fullperiod1)+5*waitbetween+1:end) = 0;
    
    % long ramp with wait inbetween (vol) at 100ms duration
    % vol1 = eField(1); vol2 = eField(2); vol3 = eField(3); vol4 = eField(4);
    % vol5 = eField(5); vol6 = eField(6); dur = eField(7)
    % numofframes = eField(8); waitbetween = eField(9); tpf = eField(10)
    elseif length(eField) == 13
        period1 = eField(1)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        period2 = eField(2)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        period3 = eField(3)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        period4 = eField(4)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        period5 = eField(5)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        period6 = eField(6)*[ones(1, 100*eField(7)/eField(10)), -ones(1, 100*eField(7)/eField(10))];
        fullperiod1 = repmat(period1, [1, eField(8)/(eField(7)/eField(10)*2)]);
        fullperiod2 = repmat(period2, [1, eField(8)/(eField(7)/eField(10)*2)]);
        fullperiod3 = repmat(period3, [1, eField(8)/(eField(7)/eField(10)*2)]);
        fullperiod4 = repmat(period4, [1, eField(8)/(eField(7)/eField(10)*2)]);
        fullperiod5 = repmat(period5, [1, eField(8)/(eField(7)/eField(10)*2)]);
        fullperiod6 = repmat(period6, [1, eField(8)/(eField(7)/eField(10)*2)]);
        waitbetween = 100*eField(9)*1000/eField(10);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:100*eField(8)) = fullperiod1;
        ef(1, 100*eField(8)+1:100*eField(8)+waitbetween) = 0;
        ef(1, 100*eField(8)+waitbetween+1:200*eField(8)+waitbetween) = fullperiod2;
        ef(1, 200*eField(8)+waitbetween+1:200*eField(8)+2*waitbetween) = 0;
        ef(1, 200*eField(8)+2*waitbetween+1:300*eField(8)+2*waitbetween) = fullperiod3;
        ef(1, 300*eField(8)+2*waitbetween+1:300*eField(8)+3*waitbetween) = 0;
        ef(1, 300*eField(8)+3*waitbetween+1:400*eField(8)+3*waitbetween) = fullperiod4;
        ef(1, 400*eField(8)+3*waitbetween+1:400*eField(8)+4*waitbetween) = 0;
        ef(1, 400*eField(8)+4*waitbetween+1:500*eField(8)+4*waitbetween) = fullperiod5;
        ef(1, 500*eField(8)+6*waitbetween+1:500*eField(8)+5*waitbetween) = 0;
        ef(1, 500*eField(8)+5*waitbetween+1:600*eField(8)+5*waitbetween) = fullperiod6;
        ef(1, 600*eField(8)+5*waitbetween+1:end) = 0;
    
    % rotating efield axis
    % wait = eField(1); minang = eField(2); maxang = eField(3); degstep = eField(4);
    % dur = eField(5); cycles = eField(6); tpf = eField(7):
    % vol = eField(8); freq = eField(9)
    elseif length(eField) == 9     
        numsteps = (eField(3)-eField(2))/eField(4)+1;
        waitframes = 100*eField(1)*1000/eField(7);
        period = eField(8)*[ones(1, 100*eField(5)/eField(7)), -ones(1, 100*eField(5)/eField(7))];
        fullperiod = repmat(period, [1, eField(6)]);
        ef = zeros(1, 100*length(Theta_cum_D));
        ef(1, 1:waitframes) = 0;
        ef(1, waitframes+1:waitframes+numsteps*length(fullperiod)) = repmat(fullperiod, [1, numsteps]);
        ef(1, waitframes+numsteps*length(fullperiod)+1:end)=0;
        
        %stairfunction for rotation
        stairsx = linspace(waitframes, 100*numsteps*(eField(6)*2*eField(5)/eField(7))+waitframes, numsteps+1);
        stairsy = (eField(2):eField(4):eField(3));
        %stairsx = [stairsx, stairsx(end)+(eField(6)*2*eField(5)/eField(7))];
        stairsy = [stairsy, stairsy(end)];
        
    end
        
   
    
%% plot the data    
f2 = figure('units','normalized','outerposition',[0 0 1 1]);
f2.Name= (['spot ' num2str(i) ' of ' num2str(K)]);
    
ax1 = subplot(6, 9, [1:3, 10:12, 19:21]);
%show fit localization image
    scatter_kde(xydata(1:10:end,1),(xydata(1:10:end,2)), 'filled', 'MarkerSize', 8); hold on
    xlim([new_mean(1)-6, new_mean(1)+6]);
    ylim([new_mean(2)-6, new_mean(2)+6]);
    title('fitted localization of all frames')
    
ax2 = subplot(6, 9, [28:30, 37:39, 46:48]);
%show assigned clusters   
    plot(xydata(class==1,1),xydata(class==1,2),'b.', 'MarkerSize', 2), hold on 
    plot(xydata(class==2,1),xydata(class==2,2),'r.', 'MarkerSize', 2)
    plot(xydata(class==3,1),xydata(class==3,2),'g.', 'MarkerSize', 2)
    plot(xydata(class==4,1),xydata(class==4,2), '.', 'MarkerEdgeColor', '#00c6ff', 'MarkerSize', 2)
    plot(xydata(class==5,1),xydata(class==5,2), '.', 'MarkerEdgeColor', '#ff7800', 'MarkerSize', 2)
    plot(xydata(class==6,1),xydata(class==6,2), '.', 'MarkerEdgeColor', '#027c1f', 'MarkerSize', 2)
    plot(new_mean(1), new_mean(2), 'k*','MarkerSize', 8);
    plot(start(:,1), start(:,2), 'kx', 'MarkerSize', 8);
    xlim([new_mean(1)-6, new_mean(1)+6]);
    ylim([new_mean(2)-6, new_mean(2)+6]);
    title('fitted localizations divided in 6 clusters');
    
    linkaxes([ax1 ax2], 'xy');

subplot(6, 9, [4:9])
    plot(Theta_D_shifted)
    xlim([0, movie_length])
    title(['Angular position theta in each frame'])

subplot(6, 9, [13:18 22:27])
%radial distribution of theta from vwmc
    hist(Theta_D_shifted, 360), hold on;
    plot(w(m_new),h(m_new), 'b*')
    plot(s, 'b')
    title('radial distribution of theta');
    xlim([0 360])
    xticks([0:60:360]);
    
    
subplot(6, 9, [31:36 40:45])
%Plot theta cumulative
    framenr = [1:length(Theta_cum_D)];
    plot(Theta_cum_D, 'black'), hold on;
    plot(framenr(class==1), Theta_cum_D(class==1),'b.');
    plot(framenr(class==2), Theta_cum_D(class==2),'r.');
    plot(framenr(class==3), Theta_cum_D(class==3),'g.');
    plot(framenr(class==4), Theta_cum_D(class==4),'.', 'MarkerEdgeColor', '#00c6ff');
    plot(framenr(class==5), Theta_cum_D(class==5),'.', 'MarkerEdgeColor', '#ff7800');
    plot(framenr(class==6), Theta_cum_D(class==6),'.', 'MarkerEdgeColor', '#027c1f');
    xlim([0, movie_length])
    title(['theta cumulative with net rotation of ' num2str(turns) ' turns with ' num2str(perc) '% fitted frames']);

    
subplot(6, 9, [49:54])
%Plot eField
yyaxis left
plot(ef, 'Linewidth', 0.1);
if length(eField) == 9
    yyaxis right
    stairs(stairsx, stairsy, 'r');
        % generate colormaps
        map_b = zeros(0, 3);
        map_r = zeros(0, 3);
        if eField(3) == 90
            vec = linspace(0,255,length(stairsx)-1);
            for l=1:length(stairsx)-1
            map_b = [map_b; 0 vec(l)/255 1];
            map_r = [map_r; 1 vec(l)/255 0];
            end
        elseif eField(3) == 180
            vec = linspace(0,255,length(stairsx)/2);
            for l = 1:length(stairsx)/2
            map_b = [map_b; 0 vec(l)/255 1];
            map_r = [map_r; 1 vec(l)/255 0];
            end
            for l = (length(stairsx)/2)-1:-1:1
            map_b = [map_b; 0 vec(l)/255 1];
            map_r = [map_r; 1 vec(l)/255 0];
            end
        end       
        hold all % now plot
        for l = 1:length(stairsx)-1
            yyaxis left
            plot([stairsx(l):stairsx(l+1)], ef(stairsx(l):stairsx(l+1)), 'Color' , map_b(l,:), 'Linewidth', 0.1, 'LineStyle', '-', 'Marker', 'none');
            yyaxis right
            stairs(stairsx(l:l+1), stairsy(l:l+1), 'Color' , map_r(l,:), 'Marker', 'none', 'LineWidth', 3, 'LineStyle', '-');
            ylim([0 180])
            yticks([0:45:180]);
        end        

end
xlim([1 100*movie_length])
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel',xt/100)
if length(eField) == 5
    title(['Field with ' num2str(eField(2)) 'ms and ' num2str(eField(3)) ' cycles']);
elseif length(eField) == 8
    title(['Field with ' num2str(eField(6)) ' frames at ' num2str(eField(2)) 'ms, ' num2str(eField(3)) 'ms, ' num2str(eField(4)) 'ms and ' num2str(eField(5)) 'ms']);
elseif length(eField) == 12
    title(['Field with ' num2str(eField(6)) ' frames at ' num2str(eField(2)) 'ms, ' num2str(eField(3)) 'ms, ' num2str(eField(4)) 'ms, ' num2str(eField(5)) 'ms, ' num2str(eField(6)) 'ms, and ' num2str(eField(7)) 'ms']);
elseif length(eField) == 7
    title(['Field with ' num2str(eField(2)) 'ms and ' num2str(eField(7)) 'x ' num2str(eField(3)) ' cycles in the orientation ' varargin{1}{1}]);
elseif length(eField) == 10
    %title(['Field with ' num2str(eField(8)) ' frames at ' num2str(eField(5)) 'V, ' num2str(eField(6)) 'V and ' num2str(eField(7)) 'V']);
    title(['Field with ' num2str(eField(8)) ' frames at ' num2str(eField(2)) 'ms, ' num2str(eField(3)) 'ms and ' num2str(eField(4)) 'ms']);
elseif length(eField) == 11
    title(['Field with ' num2str(eField(2)) ' ms and ' num2str(eField(9)) ' cycles at ' num2str(eField(3)) 'V -  ' num2str(eField(8)) 'V']);
elseif length(eField) == 13
    title(['Field with ' num2str(eField(8)) ' frames at ' num2str(eField(7)) 'ms, and at ' num2str(eField(1)) 'V, ' num2str(eField(2)) 'V, ' num2str(eField(3)) 'V, ' num2str(eField(4)) 'V, ' num2str(eField(5)) 'V and ' num2str(eField(6)) 'V']);
elseif length(eField) == 9 
    title(['Field with ' num2str(eField(8)) ' V and ' num2str(eField(9)) ' Hz rotating from ' num2str(eField(2)) '° to ' num2str(eField(3)) '° in ' num2str(eField(4)) '° steps']);
end
hold off

% savefig(['spot' num2str(i) '_analysis.fig'])
print(gcf, '-dtiff', ['spot' num2str(i) '_analysis.tiff']);