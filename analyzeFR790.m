function info = analyzeFR790(info,id,S)

% I added this comment.

% Trial Matrix (TM)
% 1: reward-assigned target (1:left, 2:right)
% 2: monkey's choice direction (1:left, 2:right)
% 3: RT from the first shape onset to the saccade
% 4: a last stim duration before the saccade
% 5: target x position
% 6: target y position
% 7: RF side (1:left, 2:right)
% 8: Tin or Tout (1:Tin, 2:Tout)
% 9: Target color configuration (Joey only, 1: Green left, 2: Green right)
% 10: Chosen target color (Joey only, 1: Green, 2: Red)
% 11-30: shape index (2-9)
% 31-50: FR during the analysis window
% 51: baseline FR
% 52: pre-saccadic FR
% 53: peri-target-onset FR (for varCE analysis)

% dbstop if error

switch id
    case 1
        Monkey = 'Eli';
        total_shape = 8;
        shape_offset = 2;
        include_trump_flag = 0;
        Red_RF_flag = 0;
        % pre_Tnd = 200ms + post_Tnd = 70ms 
        t_div = 200;
        lastAccumShape2Sac = 270;% 300;
    case 2
        Monkey = 'Joey';
        total_shape = 8;
        shape_offset = 2;
        include_trump_flag = 0;
        Red_RF_flag = 1;
        t_div = 130;
        % pre_Tnd = 130ms + post_Tnd = 60ms 
        lastAccumShape2Sac = 190;% 230; %250;
end
if isfield(info,'num_cell')
    save_name = [Monkey,'PopFr.mat'];
    if info.norm_flag
        save_name = [Monkey,'PopFrNorm.mat'];
    end
else
    save_name = [Monkey,'CellFr.mat'];
end
    
% Configuring which machine is running the code
if ismac
    matlab_dir = '~/Documents/MATLAB/';
else
    matlab_dir = '~/MATLAB/';
end

S_temp = S;

% pre-processing the trial data
if 0
    
    if isfield(info,'num_cell')
        load([Monkey,'PopBeh',num2str(lastAccumShape2Sac),'.mat'],'shapeSWOE','shapeSWOE_se');
    else
        load([Monkey,'CellBeh',num2str(lastAccumShape2Sac),'.mat'],'shapeSWOE','shapeSWOE_se');
    end
    
    shapeSWOE = [0.9,-0.9,0.7,-0.7,0.5,-0.5,0.3,-0.3]';
    
    fig = info.fig;
    
    TM = info.TM;
    FM = info.FM;
    SM = info.SM;
    KM = info.KM;
    % spM = info.spM;
    FRgain = info.FRgain;
    
    popID = info.popID;
    
    if id==2
        pick = true(size(TM,1),1);
        if Red_RF_flag % select Red in RF trials only
            pick = logical(TM(:,9)==1) & pick;
        end
        if ~include_trump_flag
            shapeMtrx = TM(:,11:30)+1;
            pick = ~logical(sum(shapeMtrx==1|shapeMtrx==2,2)) & pick;
        end
        % for Joey, analyzing only Red in RF trials (TM(:,9)==1)
        TM = TM(pick,:);
        FM = FM(pick,:);
        SM = SM(pick,:);
        KM = KM(pick,:);
        if isfield(info,'num_cell')
            popID = popID(pick);
        end
    end

    trueDeciWOE = [9999 -9999 9 -9 7 -7 5 -5 3 -3 nan];
    shapeMtrx = TM(:,11:30)+1;
    frMtrx{1} = TM(:,31:50);
    sacMtrx = SM;
    num_shown = sum(isfinite(shapeMtrx),2);
    rew_targ = TM(:,1);
    if id==2
        % 1:Green, 2:Red
        rew_color = (rew_targ~=TM(:,9))+1;
    end
    choice = TM(:,2);
    RT = TM(:,3);
    correct = TM(:,1)==TM(:,2);
    RF_side = TM(:,7);
    T_in_out = TM(:,8);
    presac_t = info.presac_t;
    postsac_t = info.postsac_t;
    shift = info.shift;
    bin_center = linspace(presac_t,postsac_t,(postsac_t-presac_t)/shift+1);
    
    if include_trump_flag
        subDeciWOE = [shapeSWOE; nan].*10;
    else
        subDeciWOE = [nan; nan; shapeSWOE; nan].*10;
    end
    
    rawShapeMtrx = shapeMtrx;
    rawShapeMtrx(isnan(rawShapeMtrx))=11;
    rawLLR{1} = trueDeciWOE(rawShapeMtrx);
    subRawLLR{1} = subDeciWOE(rawShapeMtrx);
    rawFrMtrx{1} = frMtrx{1};
    
    % Generating random order
    rng(1);
    rand_trial_order = randperm(size(TM,1));
    rng('shuffle');
    
    num_accum = num_shown;
    % Set this to 0 to include all the observed shapes (as opposed to
    % accumulated shapes)
    if 1
        t_cutoff = lastAccumShape2Sac*ones(length(num_shown),1);
        t_cutoff = t_cutoff-TM(:,4);
        pick = t_cutoff>0;
        num_accum(pick) = num_accum(pick)-1;
        while any(t_cutoff>0)
            t_cutoff = t_cutoff-250;
            pick = t_cutoff>0;
            num_accum(pick) = num_accum(pick)-1;
        end
    end

    for i = 1:size(TM,1)
        if num_accum(i) < 0
            num_accum(i) = 0;
        end
        shapeMtrx(i,num_accum(i)+1:num_shown(i)) = nan;
        frMtrx{1}(i,num_accum(i)+1:num_shown(i)) = nan;
    end
    
    shapeMtrx(isnan(shapeMtrx))=11;
    LLR{1} = trueDeciWOE(shapeMtrx);
    
    fprintf('create woe matrix with SWOE');
    
    
    subLLR{1} = subDeciWOE(shapeMtrx);
    % info.shapeSWOE = shapeSWOE;
    
    temp = [TM(:,51),frMtrx{1}];
    deltaFrMtrx{1} = diff(temp,1,2);
    temp = [TM(:,51),rawFrMtrx{1}];
    rawDeltaFrMtrx{1} = diff(temp,1,2);

    % changing the sign of WOE so that the positive WOE favors Tin
    % assuming the monkey assigned the weights symmetrically
    % to positive and negative shapes.
    switch id
        case 1
            pick = logical(TM(:,7)==1);
        case 2
            pick = logical(TM(:,9)==2);
    end
    
    LLR{1}(pick,:) = -LLR{1}(pick,:);
    cumLLR{1} = cumsum(LLR{1},2);
    
    subLLR{1}(pick,:) = -subLLR{1}(pick,:);
    subCumLLR{1} = cumsum(subLLR{1},2);
    
    rawLLR{1}(pick,:) = -rawLLR{1}(pick,:);
    rawCumLLR{1} = cumsum(rawLLR{1},2);
    
    subRawLLR{1}(pick,:) = -subRawLLR{1}(pick,:);
    subRawCumLLR{1} = cumsum(subRawLLR{1},2);
    
    LLR{2} = nan(size(LLR{1}));
    subLLR{2} = nan(size(subLLR{1}));
    cumLLR{2} = nan(size(cumLLR{1}));
    subCumLLR{2} = nan(size(subCumLLR{1}));
    frMtrx{2} = nan(size(frMtrx{1}));
    deltaFrMtrx{2} = nan(size(deltaFrMtrx{1}));
    
    rawLLR{2} = nan(size(rawLLR{1}));
    rawCumLLR{2} = nan(size(rawCumLLR{1}));
    subRawLLR{2} = nan(size(subRawLLR{1}));
    subRawCumLLR{2} = nan(size(subRawCumLLR{1}));
    rawFrMtrx{2} = nan(size(rawFrMtrx{1}));
    rawDeltaFrMtrx{2} = nan(size(rawDeltaFrMtrx{1}));
    
    % For backward analyses
    for ci = 1:20
        num_shift = 20-ci;
        % dealing with accumulated shapes
        pick = logical(num_accum==ci); 
        LLR{2}(pick,:) = fliplr(circshift(LLR{1}(pick,:),[0,num_shift]));
        cumLLR{2}(pick,:) = fliplr(circshift(cumLLR{1}(pick,:),[0,num_shift]));
        subLLR{2}(pick,:) = fliplr(circshift(subLLR{1}(pick,:),[0,num_shift]));
        subCumLLR{2}(pick,:) = fliplr(circshift(subCumLLR{1}(pick,:),[0,num_shift]));
        frMtrx{2}(pick,:) = fliplr(circshift(frMtrx{1}(pick,:),[0,num_shift]));
        deltaFrMtrx{2}(pick,:) = fliplr(circshift(deltaFrMtrx{1}(pick,:),[0,num_shift]));
        
        % dealing with shown shapes
        pick = logical(num_shown==ci);
        rawLLR{2}(pick,:) = fliplr(circshift(rawLLR{1}(pick,:),[0,num_shift]));
        rawCumLLR{2}(pick,:) = fliplr(circshift(rawCumLLR{1}(pick,:),[0,num_shift]));
        subRawLLR{2}(pick,:) = fliplr(circshift(subRawLLR{1}(pick,:),[0,num_shift]));
        subRawCumLLR{2}(pick,:) = fliplr(circshift(subRawCumLLR{1}(pick,:),[0,num_shift]));
        rawFrMtrx{2}(pick,:) = fliplr(circshift(rawFrMtrx{1}(pick,:),[0,num_shift]));
        rawDeltaFrMtrx{2}(pick,:) = fliplr(circshift(rawDeltaFrMtrx{1}(pick,:),[0,num_shift]));
    end
    
    % Removing unnecessary data
    info = rmfield(info,'TM');
    info = rmfield(info,'SM');
    info = rmfield(info,'RM');
    info = rmfield(info,'FM');
    info = rmfield(info,'KM');
    info = rmfield(info,'LM');
    % info = rmfield(info,'spM');
    clearvars S
    save_dir = [matlab_dir,'MonkeyPhys/790_sk/data/'];
    save([save_dir,save_name]);
    S = S_temp;
else
    
    % Removing unnecessary data
    % info = rmfield(info,'TM');
    % info = rmfield(info,'SM');
    % info = rmfield(info,'RM');
    % info = rmfield(info,'FM');
    % info = rmfield(info,'KM');
    % info = rmfield(info,'LM');
    save_dir = [matlab_dir,'MonkeyPhys/790_sk/data/'];
    load([save_dir,save_name]);
       
end

% causal FR is used in #11, 12, 15, 

% fig_switch:
%1: FR vs total WOE plot (include all the epochs)
%2: FR vs total WOE plot (include only one epoch at a time)
%3: Backward FR vs WOE plot (include only one epoch at a time)
%4: dFR vs dWOE plot
%5: WOE as a function of FR & t
%6: FR vs delta/cumulative WOE correlation (ANOVA) 
%7: end FR vs cumulative WOE correlation (from shape onset)
%8: Saccadic FR conditioned on N*
%9: Saccadic FR conditioned on WOE
%10: Saccadic FR conditioned on delta WOE
%11: Causal initial FR vs delta/cumulative WOE correlation
%12: Causal saccadic FR vs delta/cumulative WOE correlation (from shape onset)
%13: Causal saccadic FR conditioned on delta WOE
%14: Causal saccadic FR vs cumulative WOE correlation (time stepped)
%15: Causal saccadic FR sorted by RT
%16: FR vs total WOE partial correlation (ANCOVA)
%17: VarCE analysis
%18: Forward mixture vs non-mixture model
%19: Backward mixture vs non-mixture test
%20: Epoch-wise frequency histogram as a function of WOE & FR
%21: Variance analysis

if isfield(S,'fig_switch')
    fig_temp=S.fig_switch;
    fig_switch = zeros(1,30);
    for ai = 1:length(fig_temp);
        fig_switch(fig_temp(ai))=1;
    end
end

%% FR vs WOE plot (include all the epochs)
if fig_switch(1)
    num_fig = 1;
    uniqueWoe = unique(cumLLR{1}(isfinite(cumLLR{1})));
    pick = find(uniqueWoe>-9000 & uniqueWoe<9000);
    uniqueWoe = uniqueWoe(pick);
    for i = 1:length(uniqueWoe)
        fr_pick = frMtrx{1}(logical(cumLLR{1}==uniqueWoe(i)));
        fr_all{i} = fr_pick;
        fr_mean(i) = mean(fr_pick);
        fr_sd(i)   = std(fr_pick)/sqrt(1e3/info.window);
        fr_se(i)   = fr_sd(i)/sqrt(length(fr_pick));
    end
    
    frVec = reshape(frMtrx{1},[],1);
    cumWoeVec = reshape(cumLLR{1},[],1);
    
    pick = logical(~isnan(cumWoeVec) & cumWoeVec>-9000 & cumWoeVec<9000);
    [beta bint r rint stats]=regress(frVec(pick), [ones(sum(pick),1) cumWoeVec(pick)], 0.3173);
    slope = beta(2);
    slope_ci = (bint(2,2)-bint(2,1))/2;
    slope_se = sqrt(stats(4)./sum((cumWoeVec(pick) - mean(cumWoeVec(pick))).^2));
    p_values = stats(3);
    
    offset = beta(1);
    offset_ci = (bint(1,2)-bint(1,1))/2;
    
    maxY = nanmean(max(fr_mean))*1+20;
    x_pos = -2.5;%(uniqueWoe/10);
    
    figure(fig); clf;
    h=ploterr(uniqueWoe/10,fr_mean,[],fr_se,1,'ko','abshhy',0.01);
    set(h(1),'MarkerSize',8,'MarkerFaceColor','k')
    axis([-3 3 0 70]);
    font_size = 16;
    % text(x_pos,maxY*0.9,sprintf('Slope'),'FontSize',font_size,'FontWeight','bold');
    % text(x_pos,maxY*0.8,sprintf('%.2f\\pm%.2f(%.2f)',slope*10,slope_ci*10,p_values),'FontSize',font_size,'FontWeight','bold');
    % text(x_pos,maxY*0.7,sprintf('Offset'),'FontSize',font_size,'FontWeight','bold');
    % text(x_pos,maxY*0.6,sprintf('%.2f\\pm%.2f',offset,offset_ci),'FontSize',font_size,'FontWeight','bold');
    % xlabel('Cumulative WOE','FontSize',15,'FontWeight','bold'); 
    % ylabel('FR','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out')
    fig = fig+num_fig;
    
    % storing in the structure
    info.fr_slope = slope;
    info.fr_slope_ci = slope_ci;
    info.fr_offset = offset;
    info.fr_offset_ci = offset_ci;
    
end
%% FR vs WOE plot (include only one epoch at a time)
% Excluding the last shape to remove the saccadic contamination
if fig_switch(2)
    num_fig = 2;
    FH1 = figure(fig); clf;
    set(FH1,'position',[100 100 1400 500]);
    FH2 = figure(fig+1); clf; hold on;
    % set(FH2,'position',[50 50 1800 1200]);
    FH3 = figure(fig+2); clf; hold on;
    % set(FH3,'position',[100 100 1000 400]);
    
    uniqueNum = unique(num_accum);
    legend_str = [];
    minY = 0;
    maxY = 60;
    group = 5; % Set to 0 to plot for all unique logLR
    % color_map = colormap(jet(group));
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    fr_weak_woe = cell(1,15);
    T_in_out_pick = cell(1,15);
    fr_woe_all = cell(1,15);
    fr_axis = -20:80;
    num2analyze = 2;
    
    withoutLastShape = 0; % Set to 1 to exclude the last shape from the analysis
    clear slopes slopes_ci slopes_se p_values offsets offsets_ci
    
    
    % Linear regression
    % selectig rewarded Tin or Tout trials only
    
    
    figure(fig)
    for gi = 1:10 %uniqueNum(end)
        for i = 1:1 % Tin or Tout
            
            if withoutLastShape
                pick = logical(num_accum~=gi & ~isnan(cumLLR{1}(:,gi)) & cumLLR{1}(:,gi)>-9000 & cumLLR{1}(:,gi)<9000);
            else
                pick = logical(~isnan(cumLLR{1}(:,gi)) & cumLLR{1}(:,gi)>-9000 & cumLLR{1}(:,gi)<9000);
            end
            if sum(pick)>1
                [beta bint r rint stats]=regress(frMtrx{1}(pick,gi), [ones(sum(pick),1) cumLLR{1}(pick,gi)], 0.3173);
                slopes(gi,i) = beta(2)*10;
                slopes_ci(gi,i) = (bint(2,2)-bint(2,1))/2*10;
                slopes_se(gi,i) = sqrt(stats(4)./sum((cumLLR{1}(:,gi) - mean(cumLLR{1}(:,gi))).^2))*10;
                p_values(gi,i) = stats(3);

                offsets(gi,i) = beta(1);
                offsets_ci(gi,i) = (bint(1,2)-bint(1,1))/2;
            end
            
            if group
                min_woe_epoch = min(cumLLR{1}(:,gi))/10;
                max_woe_epoch = max(cumLLR{1}(:,gi))/10;
                frMtrxTemp = [frMtrx{1}(:,gi),cumLLR{1}(:,gi)];
                frMtrxTemp = frMtrxTemp(rand_trial_order,:);
                shuffled_num_accum = num_accum(rand_trial_order);
                if withoutLastShape
                    pick = logical(shuffled_num_accum~=gi & ~isnan(frMtrxTemp(:,end)) & frMtrxTemp(:,end)>-9000 & frMtrxTemp(:,end)<9000);
                else
                    pick = logical(~isnan(frMtrxTemp(:,end)) & frMtrxTemp(:,end)>-9000 & frMtrxTemp(:,end)<9000);
                end
                sortedFrMtrx = sortrows(frMtrxTemp(pick,:),2);
                num_per_group = floor(sum(pick)/group);
                for gi = 1:group
                    woe_mean(gi,i,gi) = mean(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,2))./10;
                    fr_woe_mean(gi,i,gi) = mean(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1));
                    fr_woe_sd(gi,i,gi) = std(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1))/sqrt(1e3/info.window);
                    fr_woe_se(gi,i,gi) = fr_woe_sd(gi,i,gi)./sqrt(num_per_group);
                    if gi<=8
                        subplot(2,5,gi);hold on
                        p=ploterr(woe_mean(gi,i,gi),fr_woe_mean(gi,i,gi),[],fr_woe_se(gi,i,gi),1,'o','abshhy',0);
                        set(p(1),'color',color_map(gi,:),'MarkerSize',12,'MarkerFaceColor',color_map(gi,:))
                        % Plot a linear regression line
                        woe_axis = min_woe_epoch:0.1:max_woe_epoch;
                        plot(squeeze(woe_mean(gi,i,:)), squeeze(woe_mean(gi,i,:))*slopes(gi,i)+offsets(gi,i),'k-');
                        % xlabel('Evidence','FontSize',24,'FontWeight','bold'); 
                        % ylabel('FR ± SD','FontSize',24,'FontWeight','bold');
                    end
                end
                xlim([-4 4])
                if gi>1
                    set(gca,'XTick',[],'YTick',[])
                end
            else
                clear fr_all fr_mean fr_sd fr_se
                uniqueWoe = unique(cumLLR{1}(isfinite(cumLLR{1}(:,gi)),gi));
                uniqueWoe = uniqueWoe(uniqueWoe>-9000 & uniqueWoe<9000);
                for gi = 1:length(uniqueWoe)
                    if withoutLastShape
                        pick = logical(cumLLR{1}(:,gi)==uniqueWoe(gi) & num_accum~=gi & cumLLR{1}(:,gi)>-9000 & cumLLR{1}(:,gi)<9000);
                    else
                        % pick = logical(cumLLR{1}(:,ei)==uniqueWoe(j) & cumLLR{1}(:,ei)>-9000 & cumLLR{1}(:,ei)<9000);
                        pick = logical(num_accum==num2analyze & cumLLR{1}(:,gi)==uniqueWoe(gi) & cumLLR{1}(:,gi)>-9000 & cumLLR{1}(:,gi)<9000);
                    end
                    
                    if sum(pick)>20
                        fr_pick = frMtrx{1}(pick,gi);
                        fr_all{gi,gi}   = fr_pick;
                        fr_mean(gi,gi)  = mean(fr_pick);
                        fr_var(gi,gi)   = var(fr_pick);
                        fr_sd(gi,gi)    = std(fr_pick)/sqrt(1e3/info.window);
                        fr_se(gi,gi)    = fr_sd(gi)/sqrt(length(fr_pick));
                        fr_woe_all{gi} = [fr_woe_all{gi};[uniqueWoe(gi)*ones(length(fr_pick),1),fr_pick]];
                    else
                        fr_pick = nan;
                        fr_all{gi,gi}   = nan;
                        fr_mean(gi,gi)  = nan;
                        fr_var(gi,gi)   = nan;
                        fr_sd(gi,gi)    = nan;
                        fr_se(gi,gi)    = nan;
                        fr_woe_all{gi} = [fr_woe_all{gi}];
                    end
                    
%                     if uniqueWoe(j)>=-10 && uniqueWoe(j)<=1
%                         fr_weak_woe{ei} = [fr_weak_woe{ei};fr_pick];
%                         T_in_out_pick{ei} = [T_in_out_pick{ei};T_in_out(pick)];
%                     end
                end
                info.uniqueWoe{gi} = uniqueWoe;
                info.fr_mean{gi} = fr_mean(gi,:);
                info.fr_sd{gi} = fr_sd(gi,:);
                
                figure(fig)
                subplot(3,4,gi); hold on
                plot(fr_mean(gi,:),fr_var(gi,:),'ko','MarkerSize',8)
                
                figure(fig+1)
                subplot(3,4,gi); hold on
                h=ploterr(uniqueWoe/10,fr_mean(gi,:),[],fr_var(gi,:)/16,1,'k.','abshhy',0); hold on
                ylim([0 60])
            end
            
            if 1
                % plot(uniqueWoe/10, uniqueWoe./10.*slopes(ei,i)+offsets(ei,i),'k--');
                if 1 && gi<=8
                    text(-2,maxY*0.9,sprintf('Slope'),'FontSize',18);
                    text(-2,maxY*0.8,sprintf('%.1f\\pm%.1f',slopes(gi),slopes_ci(gi)),'FontSize',18);
                    % text(-2,maxY*0.7,sprintf('Offset'));
                    % text(-2,maxY*0.6,sprintf('%.2f\\pm%.2f',offsets(ei),offsets_ci(ei)));
                    % title(sprintf('%d th shape (%d trials)',ei,sum(pick)),'FontSize',12);
                end
                axis([-2.5 2.5 minY maxY])
            end
            
            set(gca,'FontSize',32,'FontWeight','bold','Box','off','TickDir','out')
            % plotting just the fr means on the same figure
        end
    end
    
    if 0
        figure(fig+1);
        LLR_fixed_bin = -48:4:48;
        fr_fixed_bin = -20:4:80;
        [X,Y] = meshgrid(LLR_fixed_bin, fr_fixed_bin);
    
        for ei = 1:num2analyze
            subplot(2,4,ei)
            N = hist3(fr_woe_all{ei},{LLR_fixed_bin fr_fixed_bin });
            N_norm = N./repmat(sum(N,2),1,size(N,2));
            h=bar3(LLR_fixed_bin,N_norm);
            shading interp
            for i = 1:length(h)
                zdata = get(h(i),'ZData');
                set(h(i),'CData',zdata)
                % Add back edge color removed by interpolating shading
                set(h,'EdgeColor','k') 
            end
        end
    end
    
    figure(fig+2)
    subplot(1,2,1);hold on
    ploterr(1:8,slopes(1:8),[],slopes_ci(1:8),1,'r-.','abshhy',0.1);
    xlabel('i-th shape from the beginning','FontSize',15,'FontWeight','bold'); 
    ylabel('Slope (dFR/dWOE)','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    subplot(1,2,2);hold on
    ploterr(1:8,offsets(1:8),[],offsets_ci(1:8),1,'r-.','abshhy',0.1);
    xlabel('i-th shape from the beginning','FontSize',15,'FontWeight','bold'); 
    ylabel('Baseline FR (cumulative WOE = 0)','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')

    fig = fig+num_fig;
    
    % storing in the structure
    info.fr_slopes = slopes;
    info.fr_slopes_ci = slopes_ci;
    info.fr_offsets = offsets;
    info.fr_offsets_ci = offsets_ci;
    
end

%% Backward FR vs WOE plot (include only one epoch at a time)
if fig_switch(3)
    num_fig = 2;
    FH1 = figure(fig); clf;
    set(FH1,'position',[100 100 1000 600])
    FH2 = figure(fig+1); clf
    set(FH2,'position',[200 200 1000 600])
    uniqueNum = unique(num_accum);
    legend_str = [];
    minY = 0;
    maxY = 70;
    group = 5; % Set to 0 to plot for all unique logLR
    % color_map = colormap(jet(group));
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    clear slopes slopes_ci slopes_se p_values offsets offsets_ci
    for gi = 1:3 % uniqueNum
        
            % ANCOVA to test the difference in slopes for Tin vs Tout
            % 01/30/13
        if 0 % ei==1
            pick = logical(~isnan(cumLLR{2}(:,ei)) & cumLLR{2}(:,ei)>-9000 & cumLLR{2}(:,ei)<9000);
            [h,atab,ctab,stats] = aoctool(cumLLR{2}(pick,ei),frMtrx{2}(pick,ei),T_in_out(pick));
        end
        
        
        for i = 1:2 % Tin or Tout
            clear fr_all fr_mean fr_sd fr_se
            % Linear regression
            % selectig rewarded Tin or Tout trials only
            pick = logical(T_in_out==i & ~isnan(cumLLR{2}(:,gi)) & cumLLR{2}(:,gi)>-9000 & cumLLR{2}(:,gi)<9000);                
            
            
            if group
                % Simple plot
                min_woe_epoch = min(cumLLR{2}(:,gi))/10;
                max_woe_epoch = max(cumLLR{2}(:,gi))/10;
                frMtrxTemp = [frMtrx{2}(:,gi),cumLLR{2}(:,gi)];
                
                % Shuffle everything and use them for the following analyses
                frMtrxTemp_shuffled = frMtrxTemp(rand_trial_order,:);
                FM_shuffled = FM(rand_trial_order,:);
                T_in_out_shuffled = T_in_out(rand_trial_order);
                num_accum_shuffled = num_accum(rand_trial_order);
                
                % pick Tin or Tour trials with a valid logLR
                pick = logical(T_in_out_shuffled==i & ~isnan(frMtrxTemp_shuffled(:,end)) & frMtrxTemp_shuffled(:,end)>-9000 & frMtrxTemp_shuffled(:,end)<9000);
                
                % pick
                picked_frMtrxTemp = frMtrxTemp_shuffled(pick,:);
                picked_FM = FM_shuffled(pick,:);
                picked_num_accum = num_accum_shuffled(pick);
                
                % sort
                [sortedFrMtrx,sort_order] = sortrows(picked_frMtrxTemp,2);
                sorted_FM = picked_FM(sort_order,:);
                sorted_num_accum = picked_num_accum(sort_order);
                
                num_per_group = floor(sum(pick)/group);
                psth_back = nan(num_per_group,11);
                psth_back_woe_mean = nan(5,11,3);
                
                % Showing the analysis window
                figure(fig+1)
                subplot(2,3,3*i+1-gi);hold on
                h = fillRect([info.delay,0.5],[info.delay+250,100]);
                set(h,'LineStyle','none');
                
                for gi = 1:group
                    woe_mean(gi,i,gi) = mean(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,2))./10;
                    fr_woe_mean(gi,i,gi) = mean(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1));
                    fr_woe_sd(gi,i,gi) = std(sortedFrMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1))/sqrt(1e3/info.window);
                    fr_woe_se(gi,i,gi) = fr_woe_sd(gi,i,gi)./sqrt(num_per_group);

                    for sti = 1:num_per_group % subgroup trial index
                        ti = num_per_group*(gi-1)+sti; % trial index
                        psth_back(ti,:) = sorted_FM(ti,5*(sorted_num_accum(ti)-gi)+1:5*(sorted_num_accum(ti)-gi+2)+1);
                    end
                    
                    psth_back_woe_mean(gi,:,gi) = nanmean(psth_back);
                    psth_back_woe_sd(gi,:,gi)   = nanstd(psth_back);
                    num_entry = sum(isfinite(psth_back));
                    psth_back_woe_se(gi,:,gi)   = psth_back_woe_sd(gi,:,gi)./sqrt(num_entry);
                    figure(fig+1)
                    subplot(2,3,3*i+1-gi);hold on
                    fillTrace(0:50:500,psth_back_woe_mean(gi,:,gi),psth_back_woe_se(gi,:,gi),color_map(gi,:));
                    ylim([0 70])
                    set(gca,'FontSize',36,'FontWeight','bold','Box','OFF','TickDir','out')
                    
                    figure(fig)
                    subplot(2,3,3*i+1-gi);hold on
                    p=ploterr(woe_mean(gi,i,gi),fr_woe_mean(gi,i,gi),[],fr_woe_se(gi,i,gi),1,'o','abshhy',0);
                    set(p(1),'color',color_map(gi,:),'MarkerSize',12,'MarkerFaceColor',color_map(gi,:))
                    woe_axis = min_woe_epoch:0.1:max_woe_epoch;
                    
                    ylim([0 70])
                    set(gca,'FontSize',36,'FontWeight','bold','Box','OFF','TickDir','out')
                    
                end
                figure(fig)
                subplot(2,3,3*i+1-gi)
                plot(squeeze(woe_mean(gi,i,:)), squeeze(woe_mean(gi,i,:))*slopes(gi,i)+offsets(gi,i),'k-');
                % xlabel('Evidence','FontSize',24,'FontWeight','bold'); 
                % ylabel('FR ± SD','FontSize',24,'FontWeight','bold');
                xlim([-4 4])
                switch i
                    case 1
                        ylim([0 70])
                    case 2
                        ylim([0 70])
                end
                if gi<3
                    set(gca,'XTick',[],'YTick',[])
                end
            else
                clear fr_all fr_mean fr_sd fr_se
                uniqueWoe = unique(cumLLR{2}(isfinite(cumLLR{2}(:,gi)),gi));
                for gi = 1:length(uniqueWoe)
                    pick = logical(cumLLR{2}(:,gi)==uniqueWoe(gi) & T_in_out==i);
                    fr_pick = frMtrx{2}(pick,gi);
                    fr_all{gi} = fr_pick;
                    fr_mean(gi) = mean(fr_pick);
                    fr_sd(gi)   = std(fr_pick)/sqrt(1e3/info.window);
                    fr_se(gi)   = fr_sd(gi)/sqrt(length(fr_pick));
                end

                figure(fig)
                subplot(2,3,3*i+1-gi);hold on
                h=ploterr(uniqueWoe/10,fr_mean,[],fr_se,1,'k.','abshhy',0);
                plot(fr_mean/10, fr_mean/10*slopes(gi,i)+offsets(gi,i),'k--');
                ylim([0 70])
                set(gca,'FontSize',36,'FontWeight','bold','Box','OFF','TickDir','out')
            end
            
            if sum(pick)>1
                [beta bint r rint stats]=regress(frMtrx{2}(pick,gi), [ones(sum(pick),1) cumLLR{2}(pick,gi)], 0.3173);
                slopes(gi,i) = beta(2)*10;
                slopes_ci(gi,i) = (bint(2,2)-bint(2,1))/2*10;
                slopes_se(gi,i) = sqrt(stats(4)./sum((cumLLR{1}(:,gi) - mean(cumLLR{1}(:,gi))).^2))*10;
                p_values(gi,i) = stats(3);

                offsets(gi,i) = beta(1);
                offsets_ci(gi,i) = (bint(1,2)-bint(1,1))/2;
            end
            
            if 1
                font_size = 20;
                switch i
                    case 1
                        text(-2.5,maxY*0.2,sprintf('Slope'),'FontSize',font_size);
                        text(-2.5,maxY*0.1,sprintf('%.2f\\pm%.2f',slopes(gi,i),slopes_ci(gi,i)),'FontSize',font_size);
                    case 2
                        text(-2.5,maxY*0.9,sprintf('Slope'),'FontSize',font_size);
                        text(-2.5,maxY*0.8,sprintf('%.2f\\pm%.2f',slopes(gi,i),slopes_ci(gi,i)),'FontSize',font_size);
                end
                % text(-2.5,maxY*0.8,sprintf('%.2f\\pm%.2f(%.2f)',slopes(ei,i),slopes_ci(ei,i),p_values(ei,i)),'FontSize',font_size);
                % text(-2.5,maxY*0.7,sprintf('Offset'),'FontSize',font_size);
                % text(-2.5,maxY*0.6,sprintf('%.2f\\pm%.2f',offsets(ei,i),offsets_ci(ei,i)),'FontSize',font_size);
                title(sprintf('FR induced by up to %d th shape (%d trials)',gi,sum(pick)),'FontSize',12);
            end
            axis([-3 3 minY maxY])
            set(gca,'FontSize',36,'FontWeight','bold','Box','OFF','TickDir','out')
            % plotting just the fr means on the same figure
        end
    end
    fig = fig+num_fig;
end

%% dFR vs dWOE plot
if fig_switch(4)
    num_fig = 4;
    withoutLastShape = 0;
    FH1 = figure(fig);clf; hold on;
    set(FH1,'position',[100 100 600 600]);
    FH2 = figure(fig+1);clf; hold on;
    set(FH2,'position',[100 100 1500 800]);
    FH3 = figure(fig+2);clf; hold on;
    set(FH3,'position',[100 150 1200 800]);
    FH4 = figure(fig+3);clf; hold on;
    set(FH4,'position',[100 100 1200 800]);
    FH5 = figure(fig+4);clf; hold on;
    set(FH5,'position',[100 100 1000 400]);

    uniqueNum = unique(num_accum);
    uniqueNum = uniqueNum(uniqueNum<=12); % Restrict analysis up to 12 shapes
    woe = [-9 -7 -5 -3 3 5 7 9];
    cutoff_fr = 0;
    min_num = 1; % minimum number of shapes to be included in calculation
    for i = 1:length(woe)
        pick = logical(LLR{1}(:,min_num:end)==woe(i));
        % pick = logical(pick & frMtrx_(:,min_num:end-1)>cutoff_fr);
        % deltaFrMtrx_ = deltaFrMtrx{1}(:,min_num:end);
        delta_fr_pick = deltaFrMtrx{1}(pick)/4; % SPIKE COUNTS!!
        delta_fr_all{i} = delta_fr_pick;
        delta_fr_mean(i) = mean(delta_fr_pick);
        delta_fr_sd(i)   = std(delta_fr_pick); % /sqrt(1e3/info.window);
        delta_fr_se(i)   = delta_fr_sd(i)/sqrt(length(delta_fr_pick));
    end
    figure(fig);hold on
    ploterr(woe/10,delta_fr_mean,[],delta_fr_sd,1,'ko','abshhy',0.01);
    plot(woe/10,delta_fr_mean,'ko','MarkerFaceColor','k','MarkerSize',16);
    set(gca,'XTick',woe/10,'XTickMode','manual','XTickLabel',woe/10)
    xlim([-1 1]);
    xlabel('dWOE (logLR)','FontSize',24,'FontWeight','bold'); 
    ylabel('dFR (Hz)','FontSize',24,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out')
    hold off
    woe_axis = -10:10;
    minY = delta_fr_mean(1)*3;
    maxY = 10; % delta_fr_mean(end)*3;
    info.delta_fr_mean = delta_fr_mean;
    info.delta_fr_sd = delta_fr_sd;
    info.delta_fr_se = delta_fr_se;
    clear slopes slopes_ci slopes_se p_values offsets offsets_ci
    frMtrx_ = [TM(:,51),frMtrx{1}];
    for gi = 1:uniqueNum(end)
        for i = 1:length(woe)
            % pick = logical(LLR{1}(:,ei)==woe(i));
            % Removing the contamination of saccadic FR
            if withoutLastShape
                pick = logical(LLR{1}(:,gi)==woe(i) & num_accum~=gi & rew_targ==RF_side); % changed
            else
                pick = logical(LLR{1}(:,gi)==woe(i));
            end
            pick = logical(pick & frMtrx_(:,gi)>=cutoff_fr);
            % picking only if the previous fr was sufficiently high.
            % pick = logical(pick); % error found!! 05/06/12
            delta_fr_pick = deltaFrMtrx{1}(pick,gi)/4; % SPIKE COUNTS!!
            delta_fr_all_epoch{gi,i} = delta_fr_pick;
            delta_fr_mean_epoch(gi,i) = nanmean(delta_fr_pick);
            delta_fr_sd_epoch(gi,i)   = nanstd(delta_fr_pick); %/sqrt(1e3/info.window);
            delta_fr_se_epoch(gi,i)   = delta_fr_sd_epoch(gi,i)/sqrt(length(delta_fr_pick));
        end
        info.epoch_delta_fr_mean{gi} = delta_fr_mean_epoch(gi,:);
        info.epoch_delta_fr_sd{gi} = delta_fr_sd_epoch(gi,:);
        if withoutLastShape
            pick = logical(LLR{1}(:,gi)>-9999 & LLR{1}(:,gi)<9999 & num_accum~=gi & rew_targ==RF_side); % changed
        else
            pick = logical(LLR{1}(:,gi)>-9999 & LLR{1}(:,gi)<9999);
        end
        if sum(pick)>0
            % [beta bint r rint stats]=regress(deltaFrMtrx{1}(pick,ei), [ones(sum(pick),1) LLR{1}(pick,ei)], 0.3173);
            [beta bint r rint stats]=regress(deltaFrMtrx{1}(pick,gi)/4, [ones(sum(pick),1) LLR{1}(pick,gi)], 0.3173);
            slopes(gi)=beta(2)*10;
            slopes_ci(gi)=(bint(2,2)-bint(2,1))/2*10;
            slopes_se(gi) = sqrt(stats(4)./sum((LLR{1}(:,gi) - mean(LLR{1}(:,gi))).^2))*10;
            p_values(gi)=stats(3);
        
            offsets(gi) = beta(1);
            offsets_ci(gi) = (bint(1,2)-bint(1,1))/2;

            figure(fig+1)
            subplot(3,5,gi); hold on
            h = ploterr(woe/10,delta_fr_mean_epoch(gi,:),[],delta_fr_sd_epoch(gi,:),1,'ko','abshhy',0.01);
            set(h(1),'MarkerFaceColor','k','MarkerSize',6)
            plot(woe_axis/10, woe_axis*beta(2)+beta(1),'k--');
            % set(gca,'XTick',woe/10,'XTickMode','manual','XTickLabel',woe/10)
            if 0
                text(-0.8,maxY*0.9,sprintf('Slope'));
                text(-0.8,maxY*0.8,sprintf('%.2f\\pm%.2f(%.2f)',slopes(ei),slopes_ci(ei),p_values(ei)));
                text(-0.8,maxY*0.7,sprintf('Offset'));
                text(-0.8,maxY*0.6,sprintf('%.2f\\pm%.2f(%.2f)',offsets(ei),offsets_ci(ei)));
                title(sprintf('FR change by the %d th shape (%d trials)',ei,sum(pick)),'FontSize',12);
            end
            % xlabel('Assigned WOE','FontSize',12); 
            % ylabel('Change in FR','FontSize',12);
            set(gca,'FontSize',32,'FontWeight','bold','Box','OFF','TickDir','out')
            if gi>1
                set(gca,'XTick',[],'YTick',[])
            end
            xlim([-1 1]);
            ylim([-10 15])
            % if maxY>minY
            %     ylim([minY maxY]);
            % else
            %     ylim([maxY minY]);
            % end
        end
        hold off
    end
    
    if 0
        dfr_bin = -50:2:100;
        for ei = 1:uniqueNum(end)
            fr_all = [];
            for i = 1:length(woe)
                dfr_all = [fr_all;delta_fr_all_epoch{ei,i}];
            end
            % subplot(length(woe),length(uniqueNum),length(uniqueNum)*(i-1)+ei)
            subplot(3,4,ei)
            hist_dfr_all = histc(dfr_all,dfr_bin);
            bar(dfr_bin,hist_dfr_all);
            xlim([-50 100])
        end
    end
    
    % Now let's look at just the last accumulated shpe.
    % Does the FR change induced by the last shape during Tin correct trials still reflect WOE?
    % Or does it just go through a stereotyped pattern?
    use_Nstar = 1;
    for gi = 8 % monkey #1 (Eli) accumulated 3-8 shapes in most of the trials. 
        for i = 1:length(woe)
            if use_Nstar
                pick = logical(num_accum<=gi & LLR{2}(:,1)==woe(i) & T_in_out==1);
                last_delta_fr_pick = deltaFrMtrx{2}(pick,1);
            else
                pick = logical(num_shown<=gi & T_in_out==1 & rawLLR{2}(:,1)==woe(i));
                last_delta_fr_pick = rawDeltaFrMtrx{2}(pick,1);
            end
            last_delta_fr_all{gi,i}  = last_delta_fr_pick;
            last_delta_fr_mean(gi,i) = mean(last_delta_fr_pick);
            last_delta_fr_sd(gi,i)   = std(last_delta_fr_pick)/sqrt(1e3/info.window);
            last_delta_fr_se(gi,i)   = last_delta_fr_sd(gi,i)/sqrt(length(last_delta_fr_pick));
        end
        if use_Nstar
            pick = logical(num_accum<=gi & T_in_out==1);
            [beta bint r rint stats]=regress(deltaFrMtrx{2}(pick,1), [ones(sum(pick),1) LLR{2}(pick,1)], 0.3173);
        else
            pick = logical(num_shown<=gi & T_in_out==1);
            [beta bint r rint stats]=regress(rawDeltaFrMtrx{2}(pick,1), [ones(sum(pick),1) rawLLR{2}(pick,1)], 0.3173);
        end
        last_slopes(gi)=beta(2)*10;
        last_slopes_ci(gi)=(bint(2,2)-bint(2,1))/2*10;
        last_slopes_se(gi) = sqrt(stats(4)./sum((LLR{1}(:,gi) - mean(LLR{1}(:,gi))).^2))*10;
        last_p_values(gi)=stats(3);
        
        last_offsets(gi) = beta(1);
        last_offsets_ci(gi) = (bint(1,2)-bint(1,1))/2;

        figure(fig+3)
        subplot(2,3,gi-2); hold on
        ploterr(woe/10,last_delta_fr_mean(gi,:),[],last_delta_fr_se(gi,:),1,'r-.','abshhy',0.01);
        plot(woe_axis/10, woe_axis*beta(2)+beta(1),'g--');
        text(-0.8,maxY*0.9,sprintf('Slope'));
        text(-0.8,maxY*0.8,sprintf('%.2f\\pm%.2f(%.2f)',last_slopes(gi),last_slopes_ci(gi),last_p_values(gi)));
        text(-0.8,maxY*0.7,sprintf('Offset'));
        text(-0.8,maxY*0.6,sprintf('%.2f\\pm%.2f(%.2f)',last_offsets(gi),last_offsets_ci(gi)));
        title(sprintf('%d-shape trials (%d trials)',gi,sum(pick)),'FontSize',12);
        xlabel('Assigned WOE','FontSize',12); 
        ylabel('Change in FR','FontSize',12);
        xlim([-1 1]);
        if maxY>minY
            % ylim([minY maxY]);
        else
            % ylim([maxY minY]);
        end
        hold off
    end
    supertitle('FR change by the last shape');
    
    figure(fig+4)
    subplot(1,2,1);hold on
    ploterr(1:7,slopes(1:7),[],slopes_ci(1:7),1,'r-.','abshhy',0.1);
    xlabel('i-th shape from the beginning','FontSize',15,'FontWeight','bold'); 
    ylabel('dFR/dWOE','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    subplot(1,2,2);hold on
    ploterr(1:7,offsets(1:7),[],offsets_ci(1:7),1,'r-.','abshhy',0.1);
    xlabel('i-th shape from the beginning','FontSize',15,'FontWeight','bold'); 
    ylabel('dFR at dWOE = 0','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    fig = fig+num_fig;
    
    % storing in the structure
    info.dfr_slopes = slopes;
    info.dfr_slopes_ci = slopes_ci;
    info.dfr_offsets = offsets;
    info.dfr_offsets_ci = offsets_ci;
    
    info.last_dfr_slopes = last_slopes;
    info.last_dfr_slopes_ci = last_slopes_ci;
    info.last_dfr_offsets = last_offsets;
    info.last_dfr_offsets_ci = last_offsets_ci;
    
end

%% WOE as a function of FR & t
if fig_switch(5)
    num_fig = 1;
    fr_bin_size = 4;
    uniqueNum = unique(num_accum);
    for ui = 1:20
        if sum(num_accum<=uniqueNum(ui))/length(num_accum) > 0.9
            break
        end
    end
    maxFR = nanmean(max(frMtrx{1},[],2))*1+20; %+std(max(frMtrx{1},2)); % global max FR
    num_bin = ceil((maxFR+fr_bin_size/2)/fr_bin_size);
    FR_bin_center = 0:fr_bin_size:(fr_bin_size*(num_bin-1));
    FR_bin_edge = -(fr_bin_size/2):fr_bin_size:(fr_bin_size*num_bin);
    for gi = 1:uniqueNum(ui)
        for i = 1:num_bin
            pick = logical(frMtrx{1}(:,gi)>FR_bin_edge(i) & frMtrx{1}(:,gi)<FR_bin_edge(i+1));
            woe_pick = LLR{1}(pick,gi);
            woe_pick = woe_pick(woe_pick>-9999 & woe_pick<9999);
            woe_all{i,gi} = woe_pick;
            woe_mean(i,gi) = mean(woe_pick);
            woe_sd(i,gi)   = std(woe_pick);
        end
    end
    FH=figure(fig);clf;
    woe_mean = woe_mean(:,1:uniqueNum(ui));
    woe_mean(isnan(woe_mean))=0;
    g = woe_mean;
    [m,n] = size(woe_mean);
    h = d2gauss(m,n,1,1);
    M = imfilter(g,h,'replicate');
    % imagesc(M,'Ydata',[0 (num_bin-1)*4],[-3 3]);
    imagesc(M,'Ydata',[0 (num_bin-1)*4],[-3 3]);
    colorbar;
    axis xy;
    title('WOE = f(FR,t)','FontSize',15,'FontWeight','bold');
    xlabel('Number of shapes','FontSize',15,'FontWeight','bold'); 
    ylabel('FR','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    fig = fig+num_fig;
end

%% FR vs delta/cumulative WOE correlation (ANOVA) 
if fig_switch(6)
    
    FH1 = figure(fig);
    
    num_fig = 2;
    dw = 100; % ms
    dw_shift = 10; % ms
    t_axis = 0:dw_shift:dw_shift*500;
    figure(info.fig);clf;
    plot_count = 1;
    oneShapeFM = nan(5e5,100);
    woeVec = nan(5e5,2);
    
    % plot the RT histogram
    if 0
        figure(fig+1);hold on;
        bar(histc(TM(:,3),0:10:3000),'k')
        shape_onset = 0:25:300;
        for si = 1:length(shape_onset)
            plot([shape_onset(si) shape_onset(si)],[0 600],'r--')
        end

        set(gca,'FontSize',12,'FontWeight','bold','Box','off','TickDir','out')
        set(gca,'XTickMode','manual','XTick',0:10:300,'XTickLabel',0:100:3000);
        return
    end
    
    % Randomize the trial order (not sorted by cells any more)
    FM_shuffled = FM(rand_trial_order,:);
    rawWoeMtrx_shuffled = rawLLR{1}(rand_trial_order,:);
    cumWoeMtrx_shuffled = cumLLR{1}(rand_trial_order,:);
    num_accum_shuffled = num_accum(rand_trial_order);
    T_in_out_shuffled = T_in_out(rand_trial_order);
    RT_shuffled = RT(rand_trial_order);
    correct_shuffled = correct(rand_trial_order);
    rew_targ_shuffled = rew_targ(rand_trial_order);
    if exist('rew_color','var')
        rew_color_shuffled = rew_color(rand_trial_order);
    end
    
    if 0
        ti = 1; % trial index
        si = 0; % shape index
        for i = 1:size(FM_shuffled,1)
            for j = 1:num_accum(i)
                oneShapeFM(ti,:) = FM_shuffled(i,25*(j-1)+1:25*(j-1)+100);
                ti = ti+1;
            end
            woeVec(si+1:si+num_accum_shuffled(i),1) = rawWoeMtrx_shuffled(i,1:num_accum_shuffled(i));
            woeVec(si+1:si+num_accum_shuffled(i)-1,2) = rawWoeMtrx_shuffled(i,2:num_accum_shuffled(i));
            si = si+num_accum_shuffled(i);
        end
        woeVec = woeVec(1:si,:);
        oneShapeFM = oneShapeFM(1:si,:);
        oneShapeDeltaFM = zeros(size(oneShapeFM));
        oneShapeDeltaFM = oneShapeFM - repmat(oneShapeFM(:,1),1,100);
    end
    firstShapeFM = FM_shuffled(:,1:500);
    % Getting a difference in FR from the 1st shape onset 
    % firstShapeFM = firstShapeFM - repmat(firstShapeFM(:,1),1,500);
    
    
    % ANOVA
    for i = 1:100
        % p_all(i) = anova1(oneShapeDeltaFM(:,i),woeVec(:,1),'off');
        p_first(i) = anova1(firstShapeFM(:,i),rawWoeMtrx_shuffled(:,1),'off');
    end
    
    % Linear Regression shows a linear relation (p<0.05)
    % Eli: 240ms~, Joey 170ms~ 
    for i = 1:100
        [beta bint r rint stats] = regress(firstShapeFM(:,i),[ones(size(rawWoeMtrx_shuffled,1),1),rawWoeMtrx_shuffled(:,1)],0.05);
        p_regress(i) = stats(3);
    end

    % Simple average
    uniqueWoe = [-9 -7 -5 -3 3 5 7 9];
    figure(fig);hold on
    group = 5;
    % color_map = colormap(jet(group));
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    num_bin = size(firstShapeFM,2);
    % plot_range = [0 200; 250 450; 500 700; 750 950; 1000 1200;];
    % plot_range = [0 200; 250 450; 500 700; 750 950; 1000 1200;];
    
    min_RT = 350;
    max_RT = 600;
    num_epoch = 4;
    
    % If group is 4 or 8, FR is sorted by individual shape logLR
    % Otherwise, FR is sorted by cumulative logLR
    if group
        for gi = 1:num_epoch
            if group==8
                firstShapeTemp = [firstShapeFM,rawWoeMtrx_shuffled(:,gi)]; % already shuffled!
                for ci = 1:2
                    for gi = 1:length(uniqueWoe)
                        pick = logical(firstShapeTemp(:,end)==uniqueWoe(gi));
                        % pick = logical(firstShapeTemp(:,end)==uniqueWoe(wi) & T_in_out_shuffled==ci & RT_shuffled>min_RT & RT_shuffled<max_RT);
                        meanFirstShapeFM(ci,gi,:) = mean(firstShapeTemp(pick,1:end-1));
                        sdFirstShapeFM(ci,gi,:) = std(firstShapeTemp(pick,1:end-1));
                        seFirstShapeFM(ci,gi,:) = sdFirstShapeFM(ci,gi,:)/sqrt(sum(pick));
                        meanRT(ci,gi) = mean(RT_shuffled(pick));
                
                        subplot(2,num_epoch,num_epoch*(ci-1)+gi);hold on
                        switch id
                            case 1
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)+1) & t_axis<=t_axis(25*(gi+1)));
                                end
                            case 2
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)+1) & t_axis<=t_axis(25*(gi+1)));
                                end
                        end
                        meanFirstShapeFMtemp = squeeze(meanFirstShapeFM(ci,gi,pick));
                        % meanFirstShapeFMtemp = meanFirstShapeFMtemp - meanFirstShapeFMtemp(1);
                        p=ploterr(t_axis(pick),meanFirstShapeFMtemp,[],0*squeeze(seFirstShapeFM(ci,gi,pick)),1,'-','abshhy',0.5);
                        set(p,'color',color_map(gi,:))
                        if gi==5 || gi>=4
                            plot([meanRT(ci,gi),meanRT(ci,gi)],[0 70],'--','color',color_map(gi,:))
                        end
                        xlabel('Time from 1st shape onset','FontSize',14,'FontWeight','bold'); 
                        ylabel('FR ± SD','FontSize',18,'FontWeight','bold');
                        xlim([min(t_axis(pick)) max(t_axis(pick))])
                        switch ci
                            case 1
                                ylim([-20 70])
                            case 2
                                ylim([-20 70])
                        end
                    end
                end
            elseif group==4
                firstShapeTemp = [firstShapeFM,rawWoeMtrx_shuffled(:,gi)]; % already shuffled!
                for ci = 1:2
                    for gi = 1:length(uniqueWoe)/2
                        switch id
                            case 1
                                pick = logical((firstShapeTemp(:,end)==uniqueWoe(gi) | firstShapeTemp(:,end)==uniqueWoe(gi+1)) & rew_targ_shuffled==ci);
                            case 2
                                pick = logical((firstShapeTemp(:,end)==uniqueWoe(gi) | firstShapeTemp(:,end)==uniqueWoe(gi+1)) & rew_color_shuffled==ci);
                        end
                        % pick = logical(firstShapeTemp(:,end)==uniqueWoe(wi) & T_in_out_shuffled==ci & RT_shuffled>min_RT & RT_shuffled<max_RT);
                        meanFirstShapeFM(ci,gi,:) = mean(firstShapeTemp(pick,1:end-1));
                        sdFirstShapeFM(ci,gi,:) = std(firstShapeTemp(pick,1:end-1));
                        seFirstShapeFM(ci,gi,:) = sdFirstShapeFM(ci,gi,:)/sqrt(sum(pick));
                        meanRT(ci,gi) = mean(RT_shuffled(pick));
                
                        subplot(2,num_epoch,num_epoch*(ci-1)+gi);hold on
                        switch id
                            case 1
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)+1) & t_axis<=t_axis(25*(gi+1)));
                                end
                            case 2
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)+1) & t_axis<=t_axis(25*(gi+1)));
                                end
                        end
                        meanFirstShapeFMtemp = squeeze(meanFirstShapeFM(ci,gi,pick));
                        %%% Getting a difference in FR from the 1st shape onset %%%
                        % meanFirstShapeFMtemp = meanFirstShapeFMtemp - meanFirstShapeFMtemp(1);
                        plot(t_axis(pick),meanFirstShapeFMtemp,'color',color_map(gi,:))
                        p=ploterr(t_axis(pick)+2-gi,meanFirstShapeFMtemp,[],0*squeeze(seFirstShapeFM(ci,gi,pick)),1,'.','abshhy',0);
                        set(p,'color',color_map(gi,:))
                        if gi==5 || gi>=4
                            plot([meanRT(ci,gi),meanRT(ci,gi)],[0 70],'--','color',color_map(gi,:))
                        end
                        xlabel('Time from 1st shape onset','FontSize',14,'FontWeight','bold'); 
                        ylabel('FR ± SD','FontSize',18,'FontWeight','bold');
                        xlim([min(t_axis(pick)) max(t_axis(pick))])
                        switch ci
                            case 1
                                ylim([-5 10])
                            case 2
                                ylim([-5 10])
                        end
                    end
                end
            else
                % group by cumulative WOE
                firstShapeTemp = [firstShapeFM,cumWoeMtrx_shuffled(:,gi)]; % already shuffled!
                for i = 1:1 % Tin or Tout
                    legend_str = [];
                    % selectig rewarded Tin or Tout trials only
                    % pick = logical(T_in_out==i & rew_targ==choice);
                    % pick = logical(T_in_out==i); % changed!!
                    pick = true(size(firstShapeTemp,1),1);
                    sortedFirstShape = sortrows(firstShapeTemp(pick,:),num_bin+1);
                    pick = ~isnan(sortedFirstShape(:,end));
                    sortedFirstShape = sortedFirstShape(pick,:);
                    num_per_group = floor(sum(pick)/group);

                    for gi = 1:group
                        woe_mean(i,gi) = mean(sortedFirstShape((num_per_group*(gi-1)+1):num_per_group*gi,end));
                        meanFirstShapeDeltaFM(i,gi,:) = mean(sortedFirstShape((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1);
                        sdFirstShapeDeltaFM(i,gi,:) = std(sortedFirstShape((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1)/sqrt(1e3/info.window);
                        sdeFirstShapeDeltaFM(i,gi,:) = sdFirstShapeDeltaFM(i,gi,:)./sqrt(num_per_group);

                        subplot(2,5,gi);hold on
                        switch id
                            case 1
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)) & t_axis<=t_axis(25*(gi+1)));
                                end
                            case 2
                                if gi==1
                                    pick = logical(t_axis<=t_axis(50));
                                else
                                    pick = logical(t_axis>=t_axis(25*(gi-1)) & t_axis<=t_axis(25*(gi+1)));
                                end
                        end
                        p=ploterr(t_axis(pick),squeeze(meanFirstShapeDeltaFM(i,gi,pick)),[],squeeze(sdeFirstShapeDeltaFM(i,gi,pick)),1,'-','abshhy',0.5);
                        set(p,'color',color_map(gi,:))
                        for k = 1:2
                            strtemp = ['                    '];
                            temp = ['total WOE = ',num2str(round(10*woe_mean(i,gi))/10)];
                            strtemp(1:length(temp)) = temp;
                            legend_str = [legend_str;strtemp];
                        end
                    end
                end
                xlabel('Time from first shape onset','FontSize',24,'FontWeight','bold'); 
                ylabel('FR ± SD','FontSize',24,'FontWeight','bold');
                xlim([-50 1650])
                switch i
                    case 1
                        ylim([0 70])
                    case 2
                        ylim([0 70])
                end
                %legend(cellstr(legend_str),'Location','NorthEastOutside');
                set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            end
        end
    else
        % creating color index
        % map = colormap;

        % interval = floor((64-1)/(group-1));
        % plot_color = map(1:interval:(1+interval*(group-1)),:);
        color_map = colormap(jet(length(uniqueWoe)));
        for i = 1:length(uniqueWoe)
            pick = logical(woeVec(:,1)==uniqueWoe(i));
            plot_color = color_map(i,:);
             
            meanOneShapeDeltaFM(i,:) = mean(oneShapeDeltaFM(pick,:),1);
            sdOneShapeDeltaFM(i,:) = std(oneShapeDeltaFM(pick,:),1)/sqrt(1e3/info.window);
            sdeOneShapeDeltaFM(i,:) = sdOneShapeDeltaFM(i,:)/sqrt(sum(pick));
            subplot(2,1,1);hold on
            h=ploterr(t_axis(1:20),meanOneShapeDeltaFM(i,:),[],sdeOneShapeDeltaFM(i,:),1,'-','abshhy',0.5);
            set(h(1),'color',plot_color);
            set(h(2),'color',plot_color);
            % title('all shapes')
            xlabel('Time from shape onset (ms)','FontSize',24,'FontWeight','bold')
            set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out')

            pick = logical(rawLLR{1}(:,1)==uniqueWoe(i));
            meanFirstShapeDeltaFM(i,:) = mean(firstShapeFM(pick,:),1);
            sdFirstShapeDeltaFM(i,:) = std(firstShapeFM(pick,:),1)/sqrt(1e3/info.window);
            sdeFirstShapeDeltaFM(i,:) = sdFirstShapeDeltaFM(i,:)/sqrt(sum(pick));
            subplot(2,1,2);hold on
            h=ploterr(t_axis(1:20),meanFirstShapeDeltaFM(i,1:20),[],sdeFirstShapeDeltaFM(i,1:20),1,'-','abshhy',0.5);
            set(h(1),'color',plot_color);
            set(h(2),'color',plot_color);
            % title('1st shape only')
            xlim([0 600])
            ylim([10 40])
            xlabel('Time from shape onset (ms)','FontSize',24,'FontWeight','bold')
            set(gca,'FontSize',18,'FontWeight','bold','Box','off')
        end
    end
    set(gcf,'position',[100 100 1600 800],'Color','w')
    return
    % Correlation
    cumWoeVec = cumsum(woeVec,2);
    R = corr(cumWoeVec(:,1),oneShapeDeltaFM,'rows','pairwise');
    
    % Partial correlation
    pick = ~isnan(cumWoeVec(:,2));
    control = cumWoeVec(pick,2);
    cumWoeVec = cumWoeVec(pick,:);
    oneShapeDeltaFM = oneShapeDeltaFM(pick,:);
    PR = partialcorr(cumWoeVec,oneShapeDeltaFM,cumWoeVec(:,2),'rows','pairwise');
    PR = PR(1,:);
        
    figure(fig+1)
    plot(t_axis(1:100),R(1,:),'-ko','MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
    xlabel('Time from shape onset (ms)','FontSize',24,'FontWeight','bold')
    set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out')
    xlim([0 1e3])
    ylim([0 0.3])
    
    figure(fig+2)
    plot(t_axis(1:100),PR(1,:),'-ko','MarkerSize',12,'MarkerFaceColor','k','LineWidth',3);
    xlabel('Time from shape onset (ms)','FontSize',24,'FontWeight','bold')
    set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out')
    xlim([0 1e3])
    
    plot_count = plot_count+1;
    info.R = R;
    info.PR = PR;
    fig = fig+num_fig;
end

%% FR vs cumulative WOE correlation (from shape onset)
if fig_switch(7)
%     FH1 = figure(fig); hold on;
%     set(FH1,'position',[100 100 1600 1200])
    num_fig = 2;
    dw = 50; % ms
    dw_shift = 10; % ms
    t_axis = 0:dw_shift:dw_shift*500;
    
    % preallocate memory 
    oneShapeFM = nan(size(FM,1),51);
    
    lastStim2Sac = TM(:,4);
    lastStim2Sac = mod(lastStim2Sac - lastAccumShape2Sac,250);

    % Randomize the trial order (not sorted by cells any more)
    FM_shuffled = FM(rand_trial_order,:);
    rawLLR_shuffled{1} = rawLLR{1}(rand_trial_order,:);
    LLR_shuffled{2} = LLR{2}(rand_trial_order,:);
    cumLLR_shuffled{1} = cumLLR{1}(rand_trial_order,:);
    cumLLR_shuffled{2} = cumLLR{2}(rand_trial_order,:);
    num_accum_shuffled = num_accum(rand_trial_order);
    T_in_out_shuffled = T_in_out(rand_trial_order);
    RT_shuffled = RT(rand_trial_order);
    correct_shuffled = correct(rand_trial_order);
    rew_targ_shuffled = rew_targ(rand_trial_order);
    lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
    
    if exist('rew_color','var')
        rew_color_shuffled = rew_color(rand_trial_order);
    end
    
    num_bin = size(oneShapeFM,2);
    sortBy = 'cumLLR';
    woe = [-9 -7 -5 -3 3 5 7 9];
    woe2analyze = 1:8;
    num_epoch = 3;
    
    % create color map
    switch sortBy
        case 'cumLLR'
            group = 5;
        case 'LLR'
            group = 8;
    end
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    for gi = 1:num_epoch
        for i = 1:size(FM_shuffled,1)
            if num_accum_shuffled(i)-gi < 0
                continue
            end
            oneShapeFM(i,1:50) = FM_shuffled(i,25*(num_accum_shuffled(i)-gi)+1:25*(num_accum_shuffled(i)-gi)+50);
            switch sortBy
                case 'cumLLR'
                    oneShapeFM(i,51) = cumLLR_shuffled{2}(i,gi);
                case 'LLR'
                    oneShapeFM(i,51) = LLR_shuffled{2}(i,gi);
            end
        end
                
        % group by cumulative WOE
        subOneShapeFM = oneShapeFM; % already shuffled!
        for i = 1:2 % Tin or Tout
            % legend_str = [];
            % selectig rewarded Tin or Tout trials only
            % pick = logical(T_in_out==i & rew_targ==choice);
            pick_in_out = logical(T_in_out_shuffled==i); % changed!!
            % pick = true(size(subOneShapeFM,1),1);
            switch sortBy
                case 'cumLLR'
                    sortedFirstShapeFM = sortrows(subOneShapeFM(pick_in_out,:),51);
                    pick = ~isnan(sortedFirstShapeFM(:,end));
                    sortedFirstShapeFM = sortedFirstShapeFM(pick,:);
                    num_per_group = floor(sum(pick)/group);
                    for gi = 1:group
                        woe_mean(i,gi) = mean(sortedFirstShapeFM((num_per_group*(gi-1)+1):num_per_group*gi,end));
                        meanFirstShapeDeltaFM(i,gi,:) = mean(sortedFirstShapeFM((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1);
                        sdFirstShapeDeltaFM(i,gi,:) = std(sortedFirstShapeFM((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1)/sqrt(1e3/info.window);
                        sdeFirstShapeDeltaFM(i,gi,:) = sdFirstShapeDeltaFM(i,gi,:)./sqrt(num_per_group);

                        subplot(2,3,3*i-gi+1);hold on
                        p=ploterr(10:10:500,squeeze(meanFirstShapeDeltaFM(i,gi,:)),[],squeeze(sdeFirstShapeDeltaFM(i,gi,:)),1,'-','abshhy',0.5);
                        set(p,'color',color_map(gi,:))
                    end
                case 'LLR'
                    for gi = woe2analyze
                        pick_woe = logical(oneShapeFM(:,end)==woe(gi));
                        pick = pick_in_out & pick_woe;
                        fr_picked = oneShapeFM(pick,1:end-1);
                        LLR_picked = oneShapeFM(pick,end);

                        sac_fr_woe_mean(i,gi,:) = mean(fr_picked);
                        sac_fr_woe_sd(i,gi,:) = std(fr_picked)/sqrt(1e3/info.window);
                        sac_fr_woe_se(i,gi,:) = sac_fr_woe_sd(i,gi,:)./sqrt(sum(pick));

                        LLR_mean(i,gi,:) = nanmean(LLR_picked);
                        LLR_sd(i,gi,:) = nanstd(LLR_picked);
                        LLR_se(i,gi,:) = LLR_sd(i,gi,:)./sqrt(sum(pick));

                        lastStim2Sac_mean(i,gi) = mean(lastStim2Sac_shuffled(pick));

                        num_accum_group{gi,gi} = num_accum_shuffled(pick);
                        
                        subplot(2,3,3*i-gi+1);hold on
                        p=ploterr(10:10:500,squeeze(sac_fr_woe_mean(i,gi,:)),[],1*squeeze(sac_fr_woe_se(i,gi,:)),1,'-','abshhy',0);
                        set(p,'color',color_map(gi,:))
                        plot([lastStim2Sac_mean(i,gi) lastStim2Sac_mean(i,gi)]+250,[0 100],'--','color',color_map(gi,:))
                    end
                    % plot subjective woe
                    shape_order = [2,4,6,8,7,5,3,1];
                    SWOE = shapeSWOE(shape_order);
                    SWOE_se = shapeSWOE_se(shape_order);

                    switch id
                        case 1
                            TinSWOE = flipud(-SWOE(1:4));
                            TinSWOE_se = flipud(SWOE_se(1:4));
                            ToutSWOE = flipud(-SWOE(5:8));
                            ToutSWOE_se = flipud(SWOE_se(5:8));
                        case 2
                            TinSWOE = SWOE(5:8);
                            TinSWOE_se = SWOE_se(5:8);
                            ToutSWOE = SWOE(1:4);
                            ToutSWOE_se = SWOE_se(1:4);
                    end
                    figure(fig+1)
                    subplot(2,1,1);hold on;
                    for i = 1:4
                        h=ploterr(woe(4+i)/10,TinSWOE(i),[],TinSWOE_se(i),1,'o','abshhy',0);
                        set(h(1),'Color',color_map(4+i,:),'MarkerFaceColor',color_map(4+i,:),'MarkerSize',12)
                        set(h(2),'Color','k')
                    end
                    xlim([0 1])
                    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
                    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')

                    subplot(2,1,2);hold on;
                    for i = 1:4
                        h=ploterr(woe(i)/10,ToutSWOE(i),[],ToutSWOE_se(i),1,'o','abshhy',0);
                        set(h(1),'Color',color_map(i,:),'MarkerFaceColor',color_map(i,:),'MarkerSize',12)
                        set(h(2),'Color','k')
                    end
                    xlim([-1 0])
                    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
                    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
                    set(gcf,'position',[1400 100 400 800],'Color','w')
            end
            switch i
                case 1
                    ylim([0 70])
                case 2
                    ylim([0 70])
            end
        end
        switch gi
            case 1
                xlabel('Time from N th shape onset','FontSize',24,'FontWeight','bold');
            otherwise
                xlabel(sprintf('Time from N-%d th shape onset',gi-1),'FontSize',24,'FontWeight','bold'); 
        end
        % xlabel('Time from first shape onset','FontSize',24,'FontWeight','bold'); 
        ylabel('FR ± sem','FontSize',24,'FontWeight','bold');
        % xlim([-50 1650])
        
        %legend(cellstr(legend_str),'Location','NorthEastOutside');
        set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
        set(gcf,'position',[100 100 1200 800],'Color','w')
    end
    
    
end

%% saccadic FR conditioned on N*
if fig_switch(8)
    num_fig = 6;
    FH8=figure(fig);clf;hold on
    set(FH8,'position',[50 200 1000 500]);
    FH9=figure(fig+1);clf;hold on
    set(FH9,'position',[100 200 1000 500]);
    FH10=figure(fig+2);clf;hold on
    set(FH10,'position',[150 200 1500 500]);
    FH11=figure(fig+3);clf;hold on
    set(FH11,'position',[200 200 1500 500]);
    FH12=figure(fig+4);clf;hold on
    set(FH12,'position',[250 200 1500 500]);
    FH13=figure(fig+5);clf;hold on
    set(FH13,'position',[300 200 1500 500]);
    uniqueNum = 2:10;%unique(num_accum);
    % creating color index
    map = colormap;
    interval = floor((64-1)/(length(uniqueNum)-1));
    plot_color = map(1:interval:(1+interval*(length(uniqueNum)-1)),:);
    % make sure these values are consistent with createTM790.m
    
    
    maxFR = max(max(sacMtrx(:,11:end)));
    FR_bin = 0:5:maxFR;
    fig_ind = 0;
    legend_str = [];
    for i = 1:length(uniqueNum)
        fig_ind = fig_ind+1;
        for gi = 1:2; % Tin or Tout
            fr_m = [];
            pick = logical(num_accum==uniqueNum(i) & T_in_out == gi & correct);
            num_trial = sum(pick);
            for k = 1:size(sacMtrx,2)
                sac_fr_all{i,gi,k} = sacMtrx(pick,k);
                sac_fr_mean(i,gi,k)= mean(sacMtrx(pick,k));
                sac_fr_sd(i,gi,k)  = std(sacMtrx(pick,k))/sqrt(1e3/info.window);
                sac_fr_se(i,gi,k)  = sac_fr_sd(i,gi,k)/sqrt(sum(pick));
                if bin_center(k)>=-500
                    fr_m = [fr_m; [bin_center(k)*ones(sum(pick),1),sacMtrx(pick,k)]];
                end
                sac_fr_binned(:,k)= histc(sac_fr_all{i,gi,k},FR_bin);
                sac_fr_prob_binned(:,k) = sac_fr_binned(:,k)./sum(sac_fr_binned(:,k));
                %plot(bin_center(k)*ones(1,length(sacMtrx(pick,k))),sacMtrx(pick,k),'o','color',plot_color(i,:));
            end
            
            % Peri-saccadic mean FR
            if gi==1
                figure(fig);
            elseif gi==2
                figure(fig+1);
            end
            p=ploterr(bin_center,squeeze(sac_fr_mean(i,gi,:)),[],squeeze(sac_fr_se(i,gi,:)),1,'-','abshhy',shift/2);
            set(p,'color',plot_color(i,:))
            
            % 2D histogram
            if gi==1
                figure(fig+2);
            elseif gi==2
                figure(fig+3);
            end
            if isempty(fr_m)
                continue
            end
            subplot(2,5,fig_ind)
            hist3(fr_m,{-500:100:300 FR_bin});
            set(gcf,'renderer','opengl');
            set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
            title(['Trilas ended with ',num2str(uniqueNum(i)),' accumulated shapes (',num2str(num_trial),' trials)'],'FontSize',10,'FontWeight','bold');
            xlabel('Time from sac','FontSize',10,'FontWeight','bold');
            ylabel('Firing rate (sp/s)','FontSize',10,'FontWeight','bold');
            set(gca,'FontSize',10,'FontWeight','bold','Box','OFF');
            view(60,60)
            %bar3(sac_fr_prob_binned(:,11:end));
            
            % Colormap
            if gi==1
                figure(fig+4);
            elseif gi==2
                figure(fig+5);
            end
            subplot(2,5,fig_ind)
            imagesc(sac_fr_prob_binned(:,11:end),'Xdata',[-500 300],'Ydata',[0 maxFR]);
            title(['Trilas ended with ',num2str(uniqueNum(i)),' accumulated shapes (',num2str(num_trial),' trials)'],'FontSize',10,'FontWeight','bold');
            xlabel('Time from saccade (ms)','FontSize',10,'FontWeight','bold');
            ylabel('Firing rate (sp/s)','FontSize',10,'FontWeight','bold');
            set(gca,'FontSize',10,'FontWeight','bold','Box','OFF');
            colorbar;
            axis xy;
            if uniqueNum(i)<10
                legend_str = [legend_str;[' ',num2str(uniqueNum(i)),' shape trials']];
            else
                legend_str = [legend_str;[num2str(uniqueNum(i)),' shape trials']];
            end
        end
    end
    figure(fig)
    legend(cellstr(legend_str),'Location','NorthEastOutside');
    xlabel('Time from saccade','FontSize',15,'FontWeight','bold'); 
    ylabel('FR ± SE','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    figure(fig+1)
    legend(cellstr(legend_str),'Location','NorthEastOutside');
    xlabel('Time from saccade','FontSize',15,'FontWeight','bold'); 
    ylabel('FR ± SE','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    
    % threshold FR across k-shape trials
    for i = 1:5 % time
        figure(fig+5+i)
        for gi = 1:2 % choice
            subplot(1,2,gi)
            ploterr(uniqueNum,squeeze(sac_fr_mean(:,gi,end-13+i)),[],squeeze(sac_fr_se(:,gi,end-14+i)),1,'-','abshhy',0.5);
            xlabel('k-shape trials','FontSize',15,'FontWeight','bold'); 
            ylabel('FR ± SE','FontSize',15,'FontWeight','bold');
            % title(num2str(bin_center(end-13+i)))
            set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
        end
    end
    info.sac_fr_mean = sac_fr_mean;
    info.sac_fr_sd = sac_fr_sd;
    info.sac_fr_se = sac_fr_se;
    fig = fig+num_fig;
end

%% saccadic FR conditioned on WOE
if fig_switch(9)
    clearvars slopes slopes_ci slopes_se p_values offsets offsets_ci
    num_fig = 3;
    FH12=figure(fig);clf;hold on
    set(FH12,'position',[100 100 600 800]);
%     FH13=figure(fig+1);clf;hold on
%     set(FH13,'position',[200 100 1000 500]);
    num_trial = size(TM,1);
    group = 5;
    
    % creating color index
    % color_map = colormap(jet(group));
    
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    num_bin = (postsac_t-presac_t)/shift+1;
    % plot_range = [200 -400; -450 -650; -700 -900];
    plot_range = [200 -350; -400 -600; -650 -850];
    % plot_range = [200 -300; -350 -550; -600 -800];
    % plot_range = [200 -250; -300 -500; -550 -750];
    % plot_range = [200 -200; -250 -450; -500 -700];
    num_accum2analyze = 4;
    
    for gi = 1:3
        switch id
            case 1
                sacMtrxTemp = [sacMtrx,cumLLR{2}(:,gi)];
                popIDtemp = popID;
            case 2
                if Red_RF_flag
                    sacMtrxTemp = [sacMtrx,cumLLR{2}(:,gi)];
                    popIDtemp = popID;
                else
                    sacMtrxTemp = [sacMtrx,cumLLR{2}(:,gi)];
                    popIDtemp = popID;
                end
        end
        sacMtrxTemp = sacMtrxTemp(rand_trial_order,:);
        T_in_out_shuffled = T_in_out(rand_trial_order);
        num_accum_shuffled = num_accum(rand_trial_order);
        popID_shuffled = popIDtemp(rand_trial_order');
        correct_shuffled = correct(rand_trial_order');
        for i = 1:2 % Tin or Tout
            legend_str = [];
            % selectig rewarded Tin or Tout trials only
            % pick = logical(T_in_out==i & rew_targ==choice);
            pick = logical(T_in_out_shuffled==i);
            % pick = logical(T_in_out_shuffled==i & correct_shuffled==1);
            % pick = logical(T_in_out_shuffled==i & num_accum_shuffled>=num_accum2analyze); % changed!! 10/08/12
            
            sortedSacMtrx = sortrows(sacMtrxTemp(pick,:),num_bin+1);
            pick = logical(~isnan(sortedSacMtrx(:,end)) & sortedSacMtrx(:,end)>-9000 & sortedSacMtrx(:,end)<9000);
            sortedSacMtrx = sortedSacMtrx(pick,:);
            num_per_group = floor(sum(pick)/group);
            
            if gi<=num_accum2analyze
                for ti = 1:length(bin_center)
                    [beta bint r rint stats]=regress(sortedSacMtrx(:,ti), [ones(size(sortedSacMtrx,1),1) sortedSacMtrx(:,end)], 0.3173);
                    slopes(ti) = beta(2)*10;
                    slopes_ci(ti) = (bint(2,2)-bint(2,1))/2*10;
                    slopes_se(ti) = sqrt(stats(4)./sum((cumLLR{2}(:,gi) - mean(cumLLR{2}(:,gi))).^2))*10;
                    p_values(ti) = stats(3);

                    offsets(ti) = beta(1);
                    offsets_ci(ti) = (bint(1,2)-bint(1,1))/2;
                end
                figure(fig+1);
                subplot(3,1,gi)
                errorbar(bin_center,slopes,slopes_ci)
                xlim([-1000 200])
                figure(fig+2);
                subplot(3,1,gi)
                plot(bin_center,p_values)
                xlim([-1000 200])
            end
            
            figure(fig);
            for gi = 1:group
                woe_mean(i,gi) = mean(sortedSacMtrx((num_per_group*(gi-1)+1):num_per_group*gi,end));
                sac_fr_woe_mean(i,gi,:) = mean(sortedSacMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1);
                sac_fr_woe_sd(i,gi,:) = std(sortedSacMtrx((num_per_group*(gi-1)+1):num_per_group*gi,1:end-1),1)/sqrt(1e3/info.window);
                sac_fr_woe_se(i,gi,:) = sac_fr_woe_sd(i,gi,:)./sqrt(num_per_group);
                popIDinGroup{gi,gi} = popID_shuffled((num_per_group*(gi-1)+1):num_per_group*gi);
                num_accum_group{gi,gi} = num_accum_shuffled((num_per_group*(gi-1)+1):num_per_group*gi);
                % subplot(3,1,ei);hold on
                subplot(2,1,i);hold on
                pick = logical(bin_center<=plot_range(gi,1) & bin_center>=plot_range(gi,2));
                % pick = true(1,length(bin_center));
                p=ploterr(bin_center(pick),squeeze(sac_fr_woe_mean(i,gi,pick)),[],squeeze(sac_fr_woe_se(i,gi,pick)),1,'-','abshhy',0.5);
                set(p,'color',color_map(gi,:))
                for k = 1:2
                    strtemp = ['                    '];
                    temp = ['total WOE = ',num2str(round(10*woe_mean(i,gi))/10)];
                    strtemp(1:length(temp)) = temp;
                    legend_str = [legend_str;strtemp];
                end
            end
            xlabel('Time from saccade','FontSize',18,'FontWeight','bold'); 
            ylabel('FR ± SD','FontSize',18,'FontWeight','bold');
            xlim([-1000 200])
            switch i
                case 1
                    ylim([0 70])
                case 2
                    ylim([0 70])
            end
            %legend(cellstr(legend_str),'Location','NorthEastOutside');
            set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
        end
    end
    info.sac_fr_woe_mean = sac_fr_woe_mean;
    info.sac_fr_woe_sd = sac_fr_woe_sd;
    info.sac_fr_woe_se = sac_fr_woe_se;
    fig = fig+num_fig;
end

%% saccadic FR conditioned on delta LLR
if fig_switch(10)
    clearvars slopes slopes_ci slopes_se p_values offsets offsets_ci
    num_fig = 3;
    FH = figure(fig);
    set(FH,'position',[550 100 400 800]);
    FH1=figure(fig+1);clf;hold on
    set(FH1,'position',[100 100 450 800]);

    num_trial = size(TM,1);
    group = 5;
    woe = [-9 -7 -5 -3 3 5 7 9];
    
    % creating color index
    % color_map = colormap(jet(group));
    
    interval = floor((64-1)/(group-1));
    map = colormap;
    % color_map = map(1:interval:(1+interval*(group-1)),:);
    color_map = (colormap(jet(8)));
    
    shape_order = [2,4,6,8,7,5,3,1];
    SWOE = shapeSWOE(shape_order);
    SWOE_se = shapeSWOE_se(shape_order);
    
    switch id
        case 1
            TinSWOE = flipud(-SWOE(1:4));
            TinSWOE_se = flipud(SWOE_se(1:4));
            ToutSWOE = flipud(-SWOE(5:8));
            ToutSWOE_se = flipud(SWOE_se(5:8));
        case 2
            TinSWOE = SWOE(5:8);
            TinSWOE_se = SWOE_se(5:8);
            ToutSWOE = SWOE(1:4);
            ToutSWOE_se = SWOE_se(1:4);
    end
    figure(fig)
    subplot(2,1,1);hold on;
    for i = 1:4
        h=ploterr(woe(4+i)/10,TinSWOE(i),[],TinSWOE_se(i),1,'o','abshhy',0);
        set(h(1),'Color',color_map(4+i,:),'MarkerFaceColor',color_map(4+i,:),'MarkerSize',12)
        set(h(2),'Color','k')
    end
    xlim([0 1])
    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
    %ylim([-5 3])
    subplot(2,1,2);hold on;
    for i = 1:4
        h=ploterr(woe(i)/10,ToutSWOE(i),[],ToutSWOE_se(i),1,'o','abshhy',0);
        set(h(1),'Color',color_map(i,:),'MarkerFaceColor',color_map(i,:),'MarkerSize',12)
        set(h(2),'Color','k')
    end
    xlim([-1 0])
    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
    %ylim([-5 3])
    
    num_bin = (postsac_t-presac_t)/shift+1;
    plot_range = [200 -800; -850 -1050; -1200 -1400];
    % plot_range = [200 -450; -500 -700; -750 -950];
    % plot_range = [200 -400; -450 -650; -700 -900];
    % plot_range = [200 -350; -400 -600; -650 -850];
    % plot_range = [200 -300; -350 -550; -600 -800];
    % plot_range = [200 -250; -300 -500; -550 -750];
    % plot_range = [200 -200; -250 -450; -500 -700];
    
    for gi = 1
        
        sacMtrxTemp = [sacMtrx,LLR{2}(:,gi)];
        lastStim2Sac = TM(:,4);
        lastStim2Sac = mod(lastStim2Sac - lastAccumShape2Sac,250);
            
        % There is no point in shuffling, but still doing just in case...
        sacMtrxTemp = sacMtrxTemp(rand_trial_order,:);
        T_in_out_shuffled = T_in_out(rand_trial_order);
        num_accum_shuffled = num_accum(rand_trial_order);
        popID_shuffled = popID(rand_trial_order');
        correct_shuffled = correct(rand_trial_order');
        cumLLR_shuffled = cumLLR{2}(rand_trial_order,:);
        subCumLLR_shuffled = subCumLLR{2}(rand_trial_order,:);
        lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
        
        
        for i = 1:2 % Tin or Tout
            legend_str = [];
            switch i
                case 1
                    woe2analyze = 1:8;
                case 2
                    woe2analyze = 1:8;
            end
            
            % selectig rewarded Tin or Tout trials only
            % pick = logical(T_in_out==i & rew_targ==choice);
            pick_in_out = logical(T_in_out_shuffled==i & num_accum_shuffled>=1);
            % pick = logical(T_in_out_shuffled==i & correct_shuffled==1);
            % pick = logical(T_in_out_shuffled==i & num_accum_shuffled>=num_accum2analyze); % changed!! 10/08/12
            
            figure(fig);
            for gi = woe2analyze
                pick_woe = logical(sacMtrxTemp(:,end)==woe(gi));
                pick = pick_in_out & pick_woe;
                fr_picked = sacMtrxTemp(pick,1:end-1);
                cumLLR_picked = subCumLLR_shuffled(pick,:);
                
                sac_fr_woe_mean(i,gi,:) = mean(fr_picked);
                sac_fr_woe_sd(i,gi,:) = std(fr_picked)/sqrt(1e3/info.window);
                sac_fr_woe_se(i,gi,:) = sac_fr_woe_sd(i,gi,:)./sqrt(sum(pick));
                            
                LLR_mean(i,gi,:) = nanmean(cumLLR_picked);
                LLR_sd(i,gi,:) = nanstd(cumLLR_picked);
                LLR_se(i,gi,:) = LLR_sd(i,gi,:)./sqrt(sum(pick));
                
                lastStim2Sac_mean(i,gi) = mean(lastStim2Sac_shuffled(pick));
                num_accum_group{gi,gi} = num_accum_shuffled(pick);
            end
        end
    end
    
    for fi = 1:2
        for i = 1:2 % Tin or Tout
            color_ind = 2;
            switch i
                case 1
                    woe2analyze = 1:8;
                case 2
                    woe2analyze = 1:8;
            end
            for gi = woe2analyze
                switch fi
                    case 1
                        figure(fig+1)
                        subplot(2,1,i);hold on
                        pick = logical(bin_center<=plot_range(gi,1) & bin_center>=plot_range(gi,2));
                        p=ploterr(bin_center(pick),squeeze(sac_fr_woe_mean(i,gi,pick)),[],1*squeeze(sac_fr_woe_se(i,gi,pick)),1,'-','abshhy',0);
                        % pick = true(1,length(bin_center));
                        set(p,'color',color_map(gi,:))
                        color_ind = color_ind+1;
                        plot([-lastStim2Sac_mean(i,gi) -lastStim2Sac_mean(i,gi)]-250,[0 100],'--','color',color_map(gi,:))
                        
                        xlabel('Time from saccade','FontSize',18,'FontWeight','bold'); 
                        ylabel('FR ± SD','FontSize',18,'FontWeight','bold');
                        set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
                        xlim([-800 200])
                        switch i
                            case 1
                                ylim([20 70])
                            case 2
                                ylim([0 50])
                        end
                    case 2
                        figure(fig+2)
                        subplot(2,1,i);hold on
                        p = ploterr(1:2,flipud(squeeze(LLR_mean(i,gi,1:2))),[],flipud(squeeze(LLR_se(i,gi,1:2))),1,'-','abshhy',0);
                        set(p,'color',color_map(gi,:))
                        color_ind = color_ind+1;
                        xlim([0 3]);
                end
            end
        end
    end
    
    info.sac_fr_woe_mean = sac_fr_woe_mean;
    info.sac_fr_woe_sd = sac_fr_woe_sd;
    info.sac_fr_woe_se = sac_fr_woe_se;
    fig = fig+num_fig;
end

%% causal FR vs delta/cumulative WOE correlation (ANOVA) 
if fig_switch(11)
    
    FH1 = figure(fig);clf;
    
    num_fig = 6;
    t_axis = 0:5000;
    plot_count = 1;
    oneShapeKM = nan(5e5,100);
    woeVec = nan(5e5,2);
    order_list = [{'baseline'},{'1st'},{'2nd'},{'3rd'},{'4th'},{'5th'},...
                  {'6th'},{'7th'},{'8th'},{'9th'},{'10th'},{'11th'},{'12th'}];
    % plot the RT histogram
    if 0
        figure(fig+1);hold on;
        bar(histc(TM(:,3),0:10:3000),'k')
        shape_onset = 0:25:300;
        for si = 1:length(shape_onset)
            plot([shape_onset(si) shape_onset(si)],[0 600],'r--')
        end

        set(gca,'FontSize',12,'FontWeight','bold','Box','off','TickDir','out')
        set(gca,'XTickMode','manual','XTick',0:10:300,'XTickLabel',0:100:3000);
        return
    end
    
    Tin_rew = logical(TM(:,1)==TM(:,7));
    
    max_num_trial = max(sum(histc(LLR{1},-9:9),2));
    
    % Analyze the variability (s.d.) in the intial firing rate
    res_init_fr = nan(size(KM,1),1);
    init_fr = mean(KM(:,1:100),2)/FRgain; % mean across the first 100ms
    last_fr = TM(:,52);
    for ci = 1:popID(end)
        pick = logical(popID==ci);
        res_init_fr(pick) = init_fr(pick) - mean(init_fr(pick));
        
        init_fr_mean(ci) = mean(init_fr(pick));
        init_fr_sd(ci)   = std(init_fr(pick));
        init_fr_snr(ci) = init_fr_mean(ci)/init_fr_sd(ci);
        
        pick = logical(popID==ci & T_in_out==1);
        res_last_fr(pick) = last_fr(pick) - mean(last_fr(pick));
        last_fr_mean(ci) = mean(last_fr(pick));
        last_fr_sd(ci)   = std(last_fr(pick));
        last_fr_snr(ci) = last_fr_mean(ci)/last_fr_sd(ci);
        
    end
    init_FR_mean = mean(init_fr);
    init_FR_sd = std(res_init_fr);
    init_FR_snr = init_FR_mean/init_FR_sd;
    
    last_FR_mean = mean(last_fr(T_in_out==1));
    last_FR_sd = std(res_last_fr);
    last_FR_snr = last_FR_mean/last_FR_sd;
    
    % Randomize the trial order (not sorted by cells any more)
    KM_shuffled = KM(rand_trial_order,:);
    if 1
        woeMtrx_shuffled = LLR{1}(rand_trial_order,:);
        rawWoeMtrx_shuffled = rawLLR{1}(rand_trial_order,:);
        cumWoeMtrx_shuffled = cumLLR{1}(rand_trial_order,:);
    else
        woeMtrx_shuffled = subLLR{1}(rand_trial_order,:);
        rawWoeMtrx_shuffled = subRawLLR{1}(rand_trial_order,:);
        cumWoeMtrx_shuffled = subCumLLR{1}(rand_trial_order,:);
    end
    num_accum_shuffled = num_accum(rand_trial_order);
    T_in_out_shuffled = T_in_out(rand_trial_order);
    RT_shuffled = RT(rand_trial_order);
    correct_shuffled = correct(rand_trial_order);
    rew_targ_shuffled = rew_targ(rand_trial_order);
    Tin_rew_shuffled = Tin_rew(rand_trial_order);
    if exist('rew_color','var')
        rew_color_shuffled = rew_color(rand_trial_order);
    end
    
    if 0
        ti = 1; % trial index
        si = 0; % shape index
        for i = 1:size(KM_shuffled,1)
            for j = 1:num_accum(i)
                oneShapeKM(ti,:) = KM_shuffled(i,25*(j-1)+1:25*(j-1)+100);
                ti = ti+1;
            end
            woeVec(si+1:si+num_accum_shuffled(i),1) = rawWoeMtrx_shuffled(i,1:num_accum_shuffled(i));
            woeVec(si+1:si+num_accum_shuffled(i)-1,2) = rawWoeMtrx_shuffled(i,2:num_accum_shuffled(i));
            si = si+num_accum_shuffled(i);
        end
        woeVec = woeVec(1:si,:);
        oneShapeKM = oneShapeKM(1:si,:);
        oneShapeDeltaKM = zeros(size(oneShapeKM));
        oneShapeDeltaKM = oneShapeKM - repmat(oneShapeKM(:,1),1,100);
    end
    firstShapeKM = KM_shuffled(:,1:3000);
    % Getting a difference in FR from the 1st shape onset 
    % firstShapeKM = firstShapeKM - repmat(firstShapeKM(:,1),1,500);
    
    if 0
        % ANOVA
        for i = 1:100
            % p_all(i) = anova1(oneShapeDeltaKM(:,i),woeVec(:,1),'off');
            p_first(i) = anova1(firstShapeKM(:,i),rawWoeMtrx_shuffled(:,1),'off');
        end
    end
    
    % Linear Regression shows a linear relation (p<0.05)
    % Eli: 240ms~, Joey 170ms~ 
    if 0
        for i = 1:100
            [beta bint r rint stats] = regress(firstShapeKM(:,i),[ones(size(rawWoeMtrx_shuffled,1),1),rawWoeMtrx_shuffled(:,1)],0.05);
            p_regress(i) = stats(3);
        end
    end

    % Simple average
    uniqueWoe = [-9 -7 -5 -3 3 5 7 9];
    figure(fig);hold on
    
    num_bin = size(firstShapeKM,2);
    plot_range = [0 200; 250 450; 500 700; 750 950; 1000 1200;];
    % plot_range = [0 200; 250 450; 500 700; 750 950; 1000 1200;];
    
    t_start = 100;
    t_end = 250;
    t_gap = 20;
    sort_list = {'LLR','cumLLR'};
    
    group = 8;
    num_epoch = 10; % 10;
    
    % color_map = colormap(jet(group));
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    % If group is 4 or 8, FR is sorted by individual shape logLR
    % Otherwise, FR is sorted by cumulative logLR
    allEpochDeltaFR = cell(1,8);
    for si = 1
        sortBy = sort_list{si};
        clear meanFirstShapeKM sdFirstShapeKM seFirstShapeKM
        figure(fig);clf;
        figure(fig+1);clf;
        figure(fig+2);clf;
        figure(fig+3);clf;
        for ei = 1:num_epoch
            switch sortBy
                case 'LLR'
                    switch group
                        case 8 % sort by individual logLR
                            t_plot_pick = logical(t_axis>=t_axis(250*(ei-1)+t_div) & t_axis<=t_axis(250*ei+t_div));
                            t_pick = logical(t_axis>=t_axis(250*(ei-1)+t_div+t_start) & t_axis<=t_axis(250*(ei-1)+t_div+t_end));
                            if ei==1
                                t_pre_pick = logical(t_axis<=t_axis(t_div));
                            else
                                t_pre_pick = logical(t_axis>=t_axis(250*(ei-2)+t_div+t_start) & t_axis<=t_axis(250*(ei-2)+t_div+t_end));
                            end

                            if isempty(firstShapeKM)
                                continue
                            elseif size(firstShapeKM,1)<2
                                deltaFR = (mean(double(firstShapeKM(t_pick)),2) - mean(double(firstShapeKM(t_pre_pick)),2))/FRgain;
                            else
                                deltaFR = (mean(double(firstShapeKM(:,t_pick)),2) - mean(double(firstShapeKM(:,t_pre_pick)),2))/FRgain;
                            end

                            pick_pos = logical(woeMtrx_shuffled(:,ei)>0);
                            pick_rew = logical(rew_targ_shuffled==1);
                            
                            pick_epoch = logical(num_accum_shuffled==ei+1 | num_accum_shuffled==ei+2 | num_accum_shuffled==ei+3);
                            pick = logical(~pick_pos & pick_epoch & pick_rew);
                            % pick = logical(pick_group);
                            if sum(pick)<10
                                continue
                            end
                            % with offset
                            [beta bint r rint stats] = regress(deltaFR(pick),[ones(sum(pick),1),woeMtrx_shuffled(pick,ei)],0.05);
                            % without offset
                            % [beta bint r rint stats] = regress(deltaFR(pick),woeMtrx_shuffled(pick,ei),0.05);

                            if length(beta)>1
                                offset(ei) = beta(1);
                                offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
                                slope(ei) = beta(2)*10;
                                slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
                                p_values(ei) = stats(3);
                            else
                                offset(ei) = 0;
                                offset_ci(ei) = 0;
                                slope(ei) = beta(1)*10;
                                slope_ci(ei) = (bint(1,2)-bint(1,1))/2*10;
                                p_values(ei) = stats(3);
                            end

                            for ci = 1
                                for wi = 1:length(uniqueWoe)

                                    meanDeltaFR(wi) = nan;
                                    sdDeltaFR(wi) = nan;
                                    seDeltaFR(wi) = nan;

                                    pick_group = logical(pick_epoch & woeMtrx_shuffled(:,ei)==uniqueWoe(wi));
                                    if sum(pick_group)<=1
                                        continue
                                    end
                                    meanFirstShapeKM(ci,wi,:) = mean(double(firstShapeKM(pick_group,:)))/FRgain;
                                    sdFirstShapeKM(ci,wi,:) = std(double(firstShapeKM(pick_group,:)))/FRgain;
                                    seFirstShapeKM(ci,wi,:) = sdFirstShapeKM(ci,wi,:)/sqrt(sum(pick_group));
                                    meanRT(ci,wi) = mean(RT_shuffled(pick_group));

                                    meanDeltaFR(wi) = mean(deltaFR(pick_group));
                                    sdDeltaFR(wi) = std(deltaFR(pick_group));
                                    seDeltaFR(wi) = sdDeltaFR(wi)/sqrt(sum(pick_group));

                                    allEpochDeltaFR{wi} = [allEpochDeltaFR{wi};deltaFR(pick_group)];

                                    % subplot(1,5,wi);
                                    figure(fig);hold on;
                                    % plot FR (mean±se)
                                    t_plot = t_axis(t_plot_pick);
                                    meanFirstShapeKM_plot = squeeze(meanFirstShapeKM(ci,wi,t_plot_pick));
                                    seFirstShapeKM_plot = squeeze(seFirstShapeKM(ci,wi,t_plot_pick));
                                    fillTrace(t_plot(1:end-t_gap),meanFirstShapeKM_plot(1:end-t_gap),seFirstShapeKM_plot(1:end-t_gap),color_map(wi,:));
                                end
                            end

                            % xlabel('Time from first shape onset','FontSize',24,'FontWeight','bold'); 
                            ylabel('Firing rate (mean ± s.e.)','FontSize',24,'FontWeight','bold');
                            xlim([0 t_axis(250*(ei+1)+t_div)])
                            % xlim([min(t_axis(pick)) max(t_axis(pick))])
                            switch i
                                case 1
                                    ylim([0 70])
                                case 2
                                    ylim([0 70])
                            end


                            figure(fig+1);
                            set(gcf,'position',[100 700 900 400]);
                            subplot(2,5,ei);hold on;
                            % subplot(2,ceil(num_epoch/2),ei);hold on;

                            for wi = 1:length(uniqueWoe)
                                h = ploterr(uniqueWoe(wi)/10,meanDeltaFR(wi),[],seDeltaFR(wi),1,'ko');
                                set(h(1),'MarkerSize',12,'color','k','MarkerFaceColor','k')
                                % set(h(1),'MarkerSize',12,'color',color_map(wi,:),'MarkerFaceColor',color_map(wi,:))
                            end
                            
                            % plot([uniqueWoe(1),uniqueWoe(end)]/10,offset(ei)+slope(ei)*[uniqueWoe(1),uniqueWoe(end)]/10,'k-','LineWidth',2)
                            plot([uniqueWoe(1),uniqueWoe(4)]/10,offset(ei)+slope(ei)*[uniqueWoe(1),uniqueWoe(4)]/10,'b-','LineWidth',2)
                            % plot([uniqueWoe(5),uniqueWoe(end)]/10,offset(ei)+slope(ei)*[uniqueWoe(5),uniqueWoe(end)]/10,'r-','LineWidth',2)
                            title(sprintf('%s to %s',order_list{ei},order_list{ei+1}),'FontSize',16,'FontWeight','bold')

                            axis([-1.1 1.1 -30 30]);
                            % axis([-1.1 1.1 -50 80]);
                            maxY = 30;
                            minY = -30;
                            text(-0.9,maxY*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                            text(-0.9,maxY*0.6,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',16,'FontWeight','bold');
                            text(0.2,minY*0.3,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
                            text(0.2,minY*0.7,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',16,'FontWeight','bold');

                            figure(fig)
                            set(gca,'FontSize',20,'FontWeight','bold','Box','off','TickDir','out')
                            set(gca,'XTickMode','manual','XTick',0:250:1500);

                            figure(fig+1)
                            set(gca,'FontSize',20,'FontWeight','bold','Box','off','TickDir','out')
                            set(gca,'XTickMode','manual','XTick',[-1,0,1]);

                        case 4 % sort by individual logLR
                            for ri = 1
                                for ei = 1:length(uniqueWoe)/2
                                    switch id
                                        case 1
                                            pick = logical((rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*ei-1) | rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*ei)));
                                            % pick = logical((rawWoeMtrx_shuffled(:,ei)==uniqueWoe(wi) | rawWoeMtrx_shuffled(:,ei)==uniqueWoe(wi+1)) & rew_targ_shuffled==ri);
                                        case 2
                                            pick = logical((rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*ei-1) | rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*ei)));
                                            % pick = logical((rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*wi-1) | rawWoeMtrx_shuffled(:,ei)==uniqueWoe(2*wi)) & rew_color_shuffled==2);
                                    end
                                    % pick = logical(firstShapeTemp(:,end)==uniqueWoe(wi) & T_in_out_shuffled==ri & RT_shuffled>min_RT & RT_shuffled<max_RT);
                                    meanFirstShapeKM(ri,ei,:) = mean(double(firstShapeKM(pick,1:end)))/100;
                                    sdFirstShapeKM(ri,ei,:) = std(double(firstShapeKM(pick,1:end)))/FRgain;
                                    seFirstShapeKM(ri,ei,:) = sdFirstShapeKM(ri,ei,:)/sqrt(sum(pick));
                                    meanRT(ri,ei) = mean(RT_shuffled(pick));

                                    % figure(fig);hold on;
                                    % subplot(2,num_epoch,num_epoch*(ri-1)+ei);hold on
                                    switch id
                                        case 1
                                            if ei==1
                                                pick = logical(t_axis<=t_axis(500));
                                            else
                                                pick = logical(t_axis>=t_axis(250*(ei-1)) & t_axis<=t_axis(250*(ei+1)));
                                            end
                                        case 2
                                            if ei==1
                                                pick = logical(t_axis<=t_axis(500));
                                            else
                                                pick = logical(t_axis>=t_axis(250*(ei-1)) & t_axis<=t_axis(250*(ei+1)));
                                            end
                                    end
                                    meanFirstShapeKMtemp = squeeze(meanFirstShapeKM(ri,ei,pick));
                                    %%% Getting a difference in FR from the 1st shape onset %%%
                                    % meanFirstShapeKMtemp = meanFirstShapeKMtemp - meanFirstShapeKMtemp(1);

                                    % plot(t_axis(pick),meanFirstShapeKMtemp,'color',color_map(wi,:))
                                    % p=ploterr(t_axis(pick)+2-wi,meanFirstShapeKMtemp,[],0*squeeze(seFirstShapeKM(ri,wi,pick)),1,'.','abshhy',0);
                                    % set(p,'color',color_map(wi,:))
                                    if ei==5 || ei>=4
                                        plot([meanRT(ri,ei),meanRT(ri,ei)],[0 70],'--','color',color_map(ei,:))
                                    end
                                    xlabel('Time from 1st shape onset (ms)','FontSize',18,'FontWeight','bold'); 
                                    ylabel('Firing rate (mean±sem)','FontSize',18,'FontWeight','bold');

                                    xlim([min(t_axis(pick)) max(t_axis(pick))])
                                    switch ri
                                        case 1
                                            % ylim([-5 10])
                                        case 2
                                            % ylim([-5 10])
                                    end
                                    if ei==1
                                        figure(fig)
                                        % subplot(2,1,ri);
                                        hold on;
                                        p=ploterr(t_axis(pick)+2-ei,meanFirstShapeKMtemp,[],1*squeeze(seFirstShapeKM(ri,ei,pick)),1,'.','abshhy',0);
                                        set(p,'color',color_map(ei,:))
                                        xlabel('Time from 1st shape onset','FontSize',30,'FontWeight','bold'); 
                                        ylabel('FR ± SD','FontSize',30,'FontWeight','bold');
                                        xlim([min(t_axis(pick)) max(t_axis(pick))])
                                        set(gca,'FontSize',30,'FontWeight','bold','Box','off','TickDir','out')
                                    end
                                end
                            end
                    end

                    info.dfr_woe(ei,:) = uniqueWoe/10;
                    info.dfr_mean(ei,:) = meanDeltaFR;
                    info.dfr_sd(ei,:) = sdDeltaFR;
                    info.dfr_se(ei,:) = seDeltaFR;
                    info.dfr_slope(ei) = slope(ei);
                    info.dfr_slope_ci(ei) = slope_ci(ei);
                    info.dfr_offset(ei) = offset(ei);
                    info.dfr_offset_ci(ei) = offset_ci(ei);

        %%        
                case 'cumLLR' % sort by cumulative logLR
                    pick_rew = logical(Tin_rew_shuffled==0 | Tin_rew_shuffled==1);
                    pick = logical(isfinite(cumWoeMtrx_shuffled(:,gi)) & pick_rew);
                    
                    if 1
                        % logLR=0.3 was shown in the epoch
                        % pick = logical(pick & woeMtrx_shuffled(:,ei+1)==5);

                        % pick = logical(pick & T_in_out_shuffled==1);
                        pick = logical(pick & num_accum_shuffled>=gi);
                        
                        % pick = logical(pick & rawWoeMtrx_shuffled(:,ei+1)==3);
                    end

                    % pick_pos = logical(pick & cumWoeMtrx_shuffled(:,ei)>0);
                    [sorted_woe,sort_order] = sortrows(cumWoeMtrx_shuffled(pick,gi));
                    sortedFirstShape = firstShapeKM(pick,:);
                    sortedFirstShape = sortedFirstShape(sort_order,:);
                    % pick = ~isnan(sortedFirstShape(:,end));
                    % sortedFirstShape = sortedFirstShape(pick,:);
                    num_per_group = floor(length(sorted_woe)/group);
                    dfr_all = [];
                    
                    % Quantiles
                    for gi = 1:group
                        pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                        woe_mean(gi) = mean(sorted_woe(pick_group,end));
                        meanFirstShapeKM(gi,:) = mean(sortedFirstShape(pick_group,1:end-1),1)/FRgain;
                        sdFirstShapeKM(gi,:) = std(double(sortedFirstShape(pick_group,1:end-1)),0,1)/FRgain;
                        seFirstShapeKM(gi,:) = sdFirstShapeKM(gi,:)./sqrt(num_per_group);

                        % subplot(1,5,ei);
                        figure(fig);hold on;
                        t_gap = 20;

                        if gi==1
                            t_pre_pick = logical(t_axis<=t_axis(t_div));
                        else
                            t_pre_pick = logical(t_axis>=t_axis(250*(gi-2)+t_div+t_start) & t_axis<=t_axis(250*(gi-2)+t_div+t_end));
                        end
                        t_pick = logical(t_axis>=t_axis(250*(gi-1)+t_div+t_start) & t_axis<=t_axis(250*(gi-1)+t_div+t_end));
                        % t_pick2 = logical(t_axis>=t_axis(250*(ei)+t_div+t_start) & t_axis<=t_axis(250*(ei)+t_div+t_end));

                        dfr = (mean(sortedFirstShape(pick_group,t_pick),2)-mean(sortedFirstShape(pick_group,t_pre_pick),2))/FRgain;
                        dfr_all = [dfr_all;dfr];
                        dfr_mean(gi) = mean(dfr);
                        dfr_sd(gi) = std(dfr);
                        dfr_se(gi) = dfr_sd(gi)/sqrt(sum(pick));

                        if gi==1
                            plot_pick = logical(t_axis<=t_axis(t_div+250-t_gap));
                        else
                            plot_pick = logical(t_axis>=t_axis(250*(gi-1)+t_div) & t_axis<=t_axis(250*gi+t_div-t_gap));
                        end

                        % plot FR (mean±se)
                        fillTrace(t_axis(plot_pick),squeeze(meanFirstShapeKM(gi,plot_pick)),squeeze(seFirstShapeKM(gi,plot_pick)),color_map(gi,:));
                    end
                    
                    % xlabel('Time from first shape onset','FontSize',24,'FontWeight','bold'); 
                    ylabel('Firing rate (mean ± s.e.)','FontSize',24,'FontWeight','bold');
                    xlim([0 t_axis(250*(gi+1)+t_div)])
                    % xlim([min(t_axis(pick)) max(t_axis(pick))])
                    switch i
                        case 1
                            ylim([0 70])
                        case 2
                            ylim([0 70])
                    end

                    if sum(pick)==0
                        continue
                    end
                    
                    
                    figure(fig+1);
                    % subplot(2,ceil(num_epoch/2),ei);hold on;
                    subplot(2,5,gi);hold on;
                    maxY = 60;
                    
                    % plot for each woe
                    min_woe = min(sorted_woe);
                    max_woe = max(sorted_woe);
                    woe_bin_size = 1;
                    woe_center = min_woe:max_woe;
                    % woe_center = (min_woe+(woe_bin_size-1)/2):(max_woe-(woe_bin_size-1)/2);
                    clearvars fr_mean fr_sd fr_se
                    for gi = 1:length(woe_center)
                        pick_group = logical(sorted_woe==woe_center(gi));
                        if sum(pick_group)>1
                            fr = mean(sortedFirstShape(pick_group,t_pick),2)/FRgain;
                            fr_mean(gi) = mean(fr);
                            fr_sd(gi) = std(fr);
                            fr_se(gi) = fr_sd(gi)/sqrt(sum(pick_group));
                        else
                            fr_mean(gi) = nan;
                            fr_sd(gi) = nan;
                            fr_se(gi) = nan;
                        end
                    end
                    
                    h = ploterr(woe_center/10,fr_mean,[],fr_se,1,'-k','abshhy',0);
                    set(h(1),'LineStyle','none');
                    
                    % plot quintiles
                    for gi = 1:group
                        % plot(woe_mean(wi)/10,mean(squeeze(meanFirstShapeKM(wi,t_pick))),'o','MarkerSize',12,'color',color_map(wi,:),'MarkerFaceColor',color_map(wi,:))
                        info.fr_woe(gi,gi)  = woe_mean(gi)/10;
                        info.fr_mean(gi,gi) = mean(squeeze(meanFirstShapeKM(gi,t_pick)));
                    end

                    % linear regression
                    [beta bint r rint stats] = regress(mean(firstShapeKM(pick,t_pick),2)/FRgain,[ones(sum(pick),1),cumWoeMtrx_shuffled(pick,gi)],0.05);
                        
                        % store fitted values
                    offset(gi) = beta(1);
                    offset_ci(gi) = (bint(1,2)-bint(1,1))/2;
                    slope(gi) = beta(2)*10;
                    slope_ci(gi) = (bint(2,2)-bint(2,1))/2*10;
                    p_values(gi) = stats(3);
                    
                    plot([woe_mean(1)/10,woe_mean(end)/10],offset(gi)+slope(gi)*[woe_mean(1),woe_mean(end)]/10,'-','LineWidth',2,'color',[0 0.5 0])
                    
                        % show values on the fig
                    axis([-2.5 2.5 0 60]);
                    text(-2,maxY*0.9,sprintf('Slope'),'FontSize',12,'FontWeight','bold');
                    text(-2,maxY*0.8,sprintf('%.1f\\pm%.1f',slope(gi),slope_ci(gi)),'FontSize',12,'FontWeight','bold');
                    text(0.5,maxY*0.3,sprintf('Offset'),'FontSize',12,'FontWeight','bold');
                    text(0.5,maxY*0.2,sprintf('%.1f\\pm%.1f',offset(gi),offset_ci(gi)),'FontSize',12,'FontWeight','bold');
                    
                    % Bilinear fit
                    % pick_pos = logical(pick & (mean(firstShapeKM(:,t_pick),2)/FRgain)>10);
                    pick_pos = logical(pick & cumWoeMtrx_shuffled(:,gi)>0);
                    pick_neg = logical(pick & cumWoeMtrx_shuffled(:,gi)<0);
                    [beta_pos bint_pos r_pos rint_pos stats_pos] = regress(mean(firstShapeKM(pick_pos,t_pick),2)/FRgain,[ones(sum(pick_pos),1),cumWoeMtrx_shuffled(pick_pos,gi)],0.05);
                    [beta_neg bint_neg r_neg rint_neg stats_neg] = regress(mean(firstShapeKM(pick_neg,t_pick),2)/FRgain,[ones(sum(pick_neg),1),cumWoeMtrx_shuffled(pick_neg,gi)],0.05);
                    
                        % store fitted values
                    offset_pos(gi) = beta_pos(1);
                    offset_ci_pos(gi) = (bint_pos(1,2)-bint_pos(1,1))/2;
                    slope_pos(gi) = beta_pos(2)*10;
                    slope_ci_pos(gi) = (bint_pos(2,2)-bint_pos(2,1))/2*10;
                    p_values_pos(gi) = stats_pos(3);
                    
                    offset_neg(gi) = beta_neg(1);
                    offset_ci_neg(gi) = (bint_neg(1,2)-bint_neg(1,1))/2;
                    slope_neg(gi) = beta_neg(2)*10;
                    slope_ci_neg(gi) = (bint_neg(2,2)-bint_neg(2,1))/2*10;
                    p_values_neg(gi) = stats_neg(3);
                    
                        % slopes for positive & negative woe
                    plot([woe_mean(1)/10,0],offset_neg(gi)+slope_neg(gi)*[woe_mean(1),0]/10,'b-','LineWidth',2)
                    plot([0,woe_mean(end)/10],offset_pos(gi)+slope_pos(gi)*[0,woe_mean(end)]/10,'r-','LineWidth',2)
                    
                        % show values on the fig
                    text(-2,maxY*0.7,sprintf('Pos Slope'),'FontSize',12,'FontWeight','bold');
                    text(-2,maxY*0.6,sprintf('%.1f\\pm%.1f',slope_pos(gi),slope_ci_pos(gi)),'FontSize',12,'FontWeight','bold');
                    text(0.5,maxY*0.1,sprintf('Pos Offset'),'FontSize',12,'FontWeight','bold');
                    text(0.5,maxY*0.0,sprintf('%.1f\\pm%.1f',offset_pos(gi),offset_ci_pos(gi)),'FontSize',12,'FontWeight','bold');

                    figure(fig+2);
                    maxY = 20;
                    % subplot(2,ceil(num_epoch/2),ei);hold on;
                    subplot(2,5,gi);hold on;
                    for gi = 1:group
                        h = ploterr(woe_mean(gi)/10,dfr_mean(gi),[],dfr_se(gi),1,'ko','abshhy',0);
                        set(h(1),'MarkerSize',12,'color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:));
                    end
                    axis([-2.5 2.5 -20 20]);
                    text(-2,maxY*0.9,sprintf('delta FR'),'FontSize',16,'FontWeight','bold');
                    text(-2,maxY*0.7,sprintf('%.1f\\pm%.1f',mean(dfr_all),std(dfr_all)),'FontSize',16,'FontWeight','bold');


                    FR_all = mean(firstShapeKM(pick,t_pick),2)/FRgain;
                    cumWoe_all = cumWoeMtrx_shuffled(pick,gi);                
                    woe_conv = -9*gi:2:9*gi;
                    fr_conv = -100:100;

                    figure(fig+3)
                    subplot(2,5,gi)
                    histmat_x1 = hist2(cumWoe_all,FR_all,woe_conv,fr_conv);
                    % imagesc(woe_conv/10,fr_conv,histmat_x1);
                    imagesc(woe_conv/10,fr_conv,histmat_x1,[0 30]);
                    axis xy
                    xlabel('WOE')
                    ylabel('FR')
                    xlim([-5 5])
                    ylim([0 100])
                    hold off

                    figure(fig+4)
                    subplot(2,5,gi)
                    woe_conv = -9:2:9;
                    fr_conv = -80:80;
                    woe_all = woeMtrx_shuffled(pick,gi);
                    dfr_all = (mean(sortedFirstShape(:,t_pick),2)-mean(sortedFirstShape(:,t_pre_pick),2))/FRgain;
                    histmat_x2 = hist2(woe_all,dfr_all,woe_conv,fr_conv);
                    imagesc(woe_conv/10,fr_conv,histmat_x2);
                    axis xy
                    xlabel('delta WOE')
                    ylabel('delta FR')
                    xlim([-1 1])
                    ylim([-80 80])
                    hold off

                    figure(fig)
                    set(gcf,'position',[100 100 1000 600],'Color','w')
                    set(gca,'FontSize',30,'FontWeight','bold','Box','off','TickDir','out')
                    set(gca,'XTickMode','manual','XTick',0:250:2500);

                    figure(fig+1)
                    set(gcf,'position',[100 700 1400 500],'Color','w')
                    set(gca,'FontSize',30,'FontWeight','bold','Box','off','TickDir','out')
                    set(gca,'XTickMode','manual','XTick',[-2,0,2]);

                    figure(fig+2)
                    set(gcf,'position',[800 700 1400 500],'Color','w')
                    set(gca,'FontSize',30,'FontWeight','bold','Box','off','TickDir','out')
                    set(gca,'XTickMode','manual','XTick',[-2,0,2]);

                    figure(fig+3)
                    set(gcf,'position',[300 50 1200 400])
                    set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
                    set(gca,'XTickMode','manual','XTick',[-2,0,2]);

                    figure(fig+4)
                    set(gcf,'position',[100 100 1200 400])
                    set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
                    set(gca,'XTickMode','manual','XTick',[-2,0,2]);
                    
                    

                    info.fr_slope(gi) = slope(gi);
                    info.fr_slope_ci(gi) = slope_ci(gi);
                    info.fr_offset(gi) = offset(gi);
                    info.fr_offset_ci(gi) = offset_ci(gi);

        %%
                case 'FR'
                    pick = isfinite(cumWoeMtrx_shuffled(:,gi));

                    if 1
                        % logLR=0.3 was shown in the epoch
                        pick = logical(pick & woeMtrx_shuffled(:,gi+1)==5);
                        % pick = logical(pick & num_accum_shuffled==ei+2);
                        % pick = logical(pick & rawWoeMtrx_shuffled(:,ei+1)==3);
                    end
                    t_pick = logical(t_axis>=t_axis(250*(gi-1)+t_div+t_start) & t_axis<=t_axis(250*(gi-1)+t_div+t_end));
                    % t_pick2 = logical(t_axis>=t_axis(250*(ei)+t_div+t_start) & t_axis<=t_axis(250*(ei)+t_div+t_end));
                    % pick_pos = logical(pick & cumWoeMtrx_shuffled(:,ei)>0);

                    [sorted_fr,sort_order] = sortrows(mean(firstShapeKM(pick,t_pick)/FRgain,2));
                    sortedFirstShape = firstShapeKM(pick,:);
                    sortedFirstShape = sortedFirstShape(sort_order,:);
                    sorted_woe = cumWoeMtrx_shuffled(pick,gi);
                    sorted_woe = sorted_woe(sort_order);
                    % pick = ~isnan(sortedFirstShape(:,end));
                    % sortedFirstShape = sortedFirstShape(pick,:);
                    num_per_group = floor(length(sorted_fr)/group);
                    for gi = 1:group
                        pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                        woe_mean(gi) = mean(sorted_woe(pick_group,end));
                        fr_mean(gi) = mean(sorted_fr(pick_group));
                        meanFirstShapeKM(gi,:) = mean(sortedFirstShape(pick_group,1:end-1),1)/FRgain;
                        sdFirstShapeKM(gi,:) = std(double(sortedFirstShape(pick_group,1:end-1)),0,1)/FRgain;
                        seFirstShapeKM(gi,:) = sdFirstShapeKM(gi,:)./sqrt(num_per_group);

                        % subplot(1,5,ei);
                        figure(fig);hold on;
                        t_gap = 20;

                        t_pick = logical(t_axis>=t_axis(250*(gi-1)+t_div+t_start) & t_axis<=t_axis(250*(gi-1)+t_div+t_end));
                        t_pick2 = logical(t_axis>=t_axis(250*(gi)+t_div+t_start) & t_axis<=t_axis(250*(gi)+t_div+t_end));

                        dfr_group = (mean(sortedFirstShape(pick_group,t_pick2),2)-mean(sortedFirstShape(pick_group,t_pick),2))/FRgain;
                        dfr_mean(gi) = mean(dfr_group);
                        dfr_sd(gi) = std(dfr_group);
                        dfr_se(gi) = dfr_sd(gi)/sqrt(sum(pick));

                        if gi==1
                            plot_pick = logical(t_axis<=t_axis(t_div+250-t_gap));
                        else
                            plot_pick = logical(t_axis>=t_axis(250*(gi-1)+t_div) & t_axis<=t_axis(250*gi+t_div-t_gap));
                        end

                        % plot FR (mean±se)
                        fillTrace(t_axis(plot_pick),squeeze(meanFirstShapeKM(gi,plot_pick)),squeeze(seFirstShapeKM(gi,plot_pick)),color_map(gi,:));
                    end

                    dfr_all = (mean(sortedFirstShape(:,t_pick2),2)-mean(sortedFirstShape(:,t_pick),2))/FRgain;

                    figure(fig);
                    subplot(2,ceil(num_epoch/2),gi);hold on;
                    for gi = 1:group
                        plot(woe_mean(gi)/10,fr_mean(gi),'o','MarkerSize',12,'color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:))
                    end
                    axis([-0.5 1 0 60]);
                    % plot([woe_mean(i,1)/10,woe_mean(i,end)/10],offset(ei)+slope(ei)*[woe_mean(i,1),woe_mean(i,end)]/10,'k-')


                    figure(fig+1);
                    subplot(2,ceil(num_epoch/2),gi);hold on;
                    maxY = 20;
                    % linear regression
                    [beta bint r rint stats] = regress(dfr_all,[ones(sum(pick),1),sorted_fr],0.05);
                    % [beta bint r rint stats] = regress(mean(firstShapeKM(pick_pos,t_pick),2)/FRgain,[ones(sum(pick_pos),1),cumWoeMtrx_shuffled(pick_pos,ei)],0.05);

                    offset(gi) = beta(1);
                    offset_ci(gi) = (bint(1,2)-bint(1,1))/2;
                    slope(gi) = beta(2);
                    slope_ci(gi) = (bint(2,2)-bint(2,1))/2;
                    p_values(gi) = stats(3);

                    plot([fr_mean(1),fr_mean(end)],offset(gi)+slope(gi)*[fr_mean(1),fr_mean(end)],'k-')

                    text(40,maxY*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                    text(40,maxY*0.6,sprintf('%.1f\\pm%.1f',slope(gi),slope_ci(gi)),'FontSize',16,'FontWeight','bold');
                    text(40,maxY*0.3,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
                    text(40,maxY*0.0,sprintf('%.1f\\pm%.1f',offset(gi),offset_ci(gi)),'FontSize',16,'FontWeight','bold');

                    subplot(2,ceil(num_epoch/2),gi);hold on;
                    for gi = 1:group
                        h = ploterr(fr_mean(gi),dfr_mean(gi),[],dfr_se(gi),1,'ko','abshhy',0);
                        set(h(1),'MarkerSize',12,'color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:));
                    end
                    axis([0 80 -20 20]);
                    text(40,-maxY*0.3,sprintf('delta FR'),'FontSize',16,'FontWeight','bold');
                    text(40,-maxY*0.6,sprintf('%.1f\\pm%.1f',mean(dfr_all),std(dfr_all)),'FontSize',16,'FontWeight','bold');
            end
        end
        
        
        if 0
            figure(fig+4);
            % fitting the change of slope & offset
            init_beta = [1,1];
            opts = optimset('MaxFunEvals',10.^5,'MaxIter',10.^5);
            [beta_fit,err,exitFlag,oput,grad,hessian] = fminunc(@(beta) linearFitErr(beta,1:length(slope),slope,slope_ci),init_beta,opts);
            slope_beta = beta_fit;
            if ~all(offset_ci==0)
                [beta_fit,err,exitFlag,oput,grad,hessian] = fminunc(@(beta) linearFitErr(beta,1:length(offset),offset,offset_ci),init_beta,opts);
                offset_beta = beta_fit;
            end

            subplot(1,2,1);hold on;
            ploterr(1:length(slope)-1,slope(1:end-1),[],slope_ci(1:end-1),1,'k-');
            plot(1:length(slope),slope_beta(1)+slope_beta(2)*(1:length(slope)),'k--')
            y = get(gca,'Ylim');
            text(1,y(1)+diff(y)*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
            text(1,y(1)+diff(y)*0.8,sprintf('%.1f',slope_beta(2)),'FontSize',16,'FontWeight','bold');
            text(5,y(1)+diff(y)*0.9,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
            text(5,y(1)+diff(y)*0.8,sprintf('%.1f',slope_beta(1)),'FontSize',16,'FontWeight','bold');
            title('slope','FontSize',24,'FontWeight','bold')
            set(gca,'FontSize',24,'FontWeight','bold')
            set(gcf,'Color','w','position',[1100 100 800 300])

            if ~all(offset_ci==0)
                subplot(1,2,2);hold on;
                ploterr(1:length(offset)-1,offset(1:end-1),[],offset_ci(1:end-1),1,'k-');
                plot(1:length(offset),offset_beta(1)+offset_beta(2)*(1:length(offset)),'k--')
                y = get(gca,'Ylim');
                text(1,y(1)+diff(y)*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                text(1,y(1)+diff(y)*0.8,sprintf('%.1f',offset_beta(2)),'FontSize',16,'FontWeight','bold');
                text(5,y(1)+diff(y)*0.9,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
                text(5,y(1)+diff(y)*0.8,sprintf('%.1f',offset_beta(1)),'FontSize',16,'FontWeight','bold');
                title('offset','FontSize',24,'FontWeight','bold')
                set(gca,'FontSize',24,'FontWeight','bold')
            end

        end
        %%
        if strcmp(sortBy,'LLR')
            allFR = [];
            allWOE = [];
            for gi = 1:length(uniqueWoe)
                allEpochFRmean(gi) = nanmean(allEpochDeltaFR{gi});
                allEpochFRsd(gi) = nanstd(allEpochDeltaFR{gi});
                allEpochFRse(gi) = allEpochFRsd(gi)/sqrt(length(allEpochDeltaFR{gi}));
                allFR = [allFR;allEpochDeltaFR{gi}];
                allWOE = [allWOE;uniqueWoe(gi)*ones(size(allEpochDeltaFR{gi}))];
            end

            figure(fig+4);clf;hold on

            figure(fig+5);clf;hold on

            pick_pos = logical(allWOE>0);

            plot_size = 24;

            for ri = 1:1
                switch ri
                    case 1
                        [beta bint r rint stats] = regress(allFR,[ones(length(allWOE),1),allWOE],0.05);
                        fit_range = 'ALL';
                    case 2
                        [beta bint r rint stats] = regress(allFR(pick_pos),[ones(sum(pick_pos),1),allWOE(pick_pos)],0.05);
                        pick = logical(uniqueWoe>0);
                        figure(fig+4)
                        % plot([0,uniqueWoe(pick)]/10,[0,uniqueWoe(pick)]/10*beta(1)*10,'k--','LineWidth',2)
                        plot(uniqueWoe(pick)/10,uniqueWoe(pick)/10*beta(2)*10+beta(1),'ks','Markersize',plot_size,'MarkerFaceColor','r')
                        figure(fig+5)
                        % plot([uniqueWoe(~pick) 0]/10,[uniqueWoe(~pick) 0]/10*-beta(1)*10,'k--','LineWidth',2)
                        % plot(uniqueWoe(~pick)/10,uniqueWoe(~pick)/10*-beta(2)*10+beta(1),'ks','Markersize',plot_size,'MarkerFaceColor','r')
                        fit_range = 'POS';
                        info.fr_pos_woe_slope = beta(1)*10;
                    case 3
                        [beta bint r rint stats] = regress(allFR(~pick_pos),[ones(sum(~pick_pos),1),allWOE(~pick_pos)],0.05);
                        pick = logical(uniqueWoe<0);
                        figure(fig+4)
                        % plot([uniqueWoe(pick) 0]/10,[uniqueWoe(pick) 0]/10*beta(1)*10,'k--','LineWidth',2)
                        % plot(uniqueWoe(pick)/10,uniqueWoe(pick)/10*beta(2)*10+beta(1),'ks','Markersize',plot_size,'MarkerFaceColor','r')
                        figure(fig+5)
                        % plot([0,uniqueWoe(~pick)]/10,[0,uniqueWoe(~pick)]/10*-beta(1)*10,'k--','LineWidth',2)
                        plot(uniqueWoe(~pick)/10,uniqueWoe(~pick)/10*-beta(2)*10+beta(1),'ks','Markersize',plot_size,'MarkerFaceColor','r')
                        fit_range = 'NEG';
                        info.fr_neg_woe_slope = beta(1)*10;
                end
                if length(beta)>1
                    offset = beta(1);
                    offset_ci = (bint(1,2)-bint(1,1))/2;
                    slope = beta(2)*10;
                    slope_ci = (bint(2,2)-bint(2,1))/2*10;
                    p_values = stats(3);
                else
                    offset = 0;
                    offset_ci = 0;
                    slope = beta(1)*10;
                    slope_ci = (bint(1,2)-bint(1,1))/2*10;
                    p_values = stats(3);
                end

                y = get(gca,'Ylim');
                if 1
                    text(-0.1,y(1)+diff(y)*(0.9-(0.2*(ri-1))),sprintf([fit_range,' Slope']),'FontSize',16,'FontWeight','bold');
                    text(0.6,y(1)+diff(y)*(0.9-(0.2*(ri-1))),sprintf('%.1f\\pm%.1f',slope,slope_ci),'FontSize',16,'FontWeight','bold');
                end

                % storing in the structure
                info.dfr_slope_all(ri)    = slope;
                info.dfr_slope_ci_all(ri) = slope_ci;
            end
            % title('slope','FontSize',18,'FontWeight','bold')
            figure(fig+4)
            h = ploterr(uniqueWoe/10,allEpochFRmean,[],allEpochFRsd,1,'ko','abshhy',0);
            set(h(1),'MarkerSize',18,'MarkerFaceColor','k');
            xlim([-1.1 1.1])
            xlabel('single shape logLR','FontSize',24,'FontWeight','bold')
            ylabel('change in FR','FontSize',24,'FontWeight','bold')
            set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
            set(gcf,'Color','w','position',[1100 100 400 400])

            figure(fig+5)
            h = ploterr(uniqueWoe/10,fliplr(allEpochFRmean),[],fliplr(allEpochFRsd),1,'ko','abshhy',0);
            set(h(1),'MarkerSize',18,'MarkerFaceColor','k');
            xlim([-1.1 1.1])
            xlabel('single shape logLR','FontSize',24,'FontWeight','bold')
            ylabel('change in FR','FontSize',24,'FontWeight','bold')
            set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
            set(gcf,'Color','w','position',[1100 100 400 400])

            % storing in the structure
            info.dfr_mean_all = allEpochFRmean;
            info.dfr_sd_all = allEpochFRsd;
            info.dfr_se_all = allEpochFRse;

        end
        %%
        if strcmp(sortBy,'cumLLR')

            figure(fig+5)
            subplot(2,3,1)
            ploterr(1:num_epoch,slope,[],slope_ci,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            ylim([-5 15])
            subplot(2,3,2)
            ploterr(1:num_epoch,slope_pos,[],slope_ci_pos,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            ylim([-5 15])
            subplot(2,3,3)
            ploterr(1:num_epoch,slope_neg,[],slope_ci_neg,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            ylim([-5 15])
            subplot(2,3,4)
            ploterr(1:num_epoch,offset,[],offset_ci,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            subplot(2,3,5)
            ploterr(1:num_epoch,offset_pos,[],offset_ci_pos,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')
            subplot(2,3,6)
            ploterr(1:num_epoch,offset_neg,[],offset_ci_neg,1,'-ko','abshhy',0)
            set(gca,'FontSize',14,'FontWeight','bold','Box','off','TickDir','out')

            figure(fig+5)
            set(gcf,'position',[400 100 1000 600])
            
            
        end
        
    end
    fig = fig+num_fig;
end

%% causal FR vs cumulative WOE correlation (backward from saccade)
if fig_switch(12)

    num_fig = 3;
    t_axis = 0:5000;
    
    % preallocate memory 
    oneShapeKM = nan(size(KM,1),51);
    
    lastStim2Sac = TM(:,4);
    pick = logical(lastStim2Sac<lastAccumShape2Sac);
    lastStim2Sac(pick) = lastStim2Sac(pick)+250;
    
    % lastStim2Sac = mod(lastStim2Sac - lastAccumShape2Sac,250);

    % Randomize the trial order (not sorted by cells any more)
    KM_shuffled = KM(rand_trial_order,:);
    rawLLR_shuffled{1} = rawLLR{1}(rand_trial_order,:);
    LLR_shuffled = LLR{2}(rand_trial_order,:);
    cumLLR_shuffled = cumLLR{2}(rand_trial_order,:);
    num_accum_shuffled = num_accum(rand_trial_order);
    T_in_out_shuffled = T_in_out(rand_trial_order);
    RT_shuffled = RT(rand_trial_order);
    correct_shuffled = correct(rand_trial_order);
    rew_targ_shuffled = rew_targ(rand_trial_order);
    subCumLLR_shuffled = subCumLLR{2}(rand_trial_order,:);
    lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
    popID_shuffled = popID(rand_trial_order);
    if exist('rew_color','var')
        rew_color_shuffled = rew_color(rand_trial_order);
    end

    num_bin = size(oneShapeKM,2);
    sortBy = 'cumLLR';
    woe = [-9 -7 -5 -3 3 5 7 9];
    woe2analyze = 1:8;
    
    
    % create color map
    switch sortBy
        case 'cumLLR'
            num_epoch = 3;
            t_offset = 500;
            group = 5;
        case 'LLR'
            num_epoch = 3;
            t_offset = 500;
            group = 8;
    end
    
    color_group = 5;
    interval = floor((64-1)/(color_group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(color_group-1)),:);
    
    t_start = 100;
    t_end = 250;
    t_gap = 20;
    t_axis = 0:5000;
    
    
    FR_all_epoch = cell(1,2); 
    cumLLR_all_epoch = cell(1,2);
    epoch_label = cell(1,2);
    
% for ni = 4:7;
    figure(fig)
    % plot saccade time range
    t_sac_mean = mean(lastStim2Sac_shuffled) + t_offset;
    t_sac_sd = std(lastStim2Sac_shuffled);
    fill(t_sac_mean+[-t_sac_sd,t_sac_sd,t_sac_sd,-t_sac_sd], [0.2 0.2 100 100],[0.9 0.9 1],'EdgeColor','w')
    for gi = 1:num_epoch
        t_pick = logical(t_axis>=t_axis(t_div+1+t_start) & t_axis<=t_axis(t_div+t_end));
        switch id
            case 1
                t_pre_pick = logical(t_axis>=t_axis(t_div+1-150) & t_axis<=t_axis(t_div));
            case 2
                t_pre_pick = logical(t_axis<=t_axis(t_div));
        end

        for i = 1:size(KM_shuffled,1)
            if num_accum_shuffled(i)-gi < 0
                continue
            end
            oneShapeKM(i,1:500) = KM_shuffled(i,250*(num_accum_shuffled(i)-gi)+1:250*(num_accum_shuffled(i)-gi)+500);
        end
        
        if 0
            % ANCOVA to compare the FR convergence for Tin vs Tout
            [h,atab,ctab,stats] = aoctool(cumLLR_shuffled(:,ei)/10,mean(double(oneShapeKM(:,t_pick))/FRgain,2),T_in_out_shuffled,[],'cumLLR','FR','T_in_out');
        end

        for i = 2:-1:1 % Tin or Tout
            clearvars epochGroupFR
            pick_num_accum = logical(num_accum_shuffled>0);
            % pick_num_accum = logical(num_accum_shuffled<=4);
            pick_in_out = logical(T_in_out_shuffled==i); % changed!!
            pick_sign = logical((-1)^i*LLR_shuffled(:,gi)<0); % choose positive for Tin and negative for Tout
            switch sortBy
                case 'cumLLR'
                    fig_h_size = 800;
                    t_plot_pick = logical(t_axis>=t_axis(250*(3-gi)+t_div) & t_axis<=t_axis(250*(4-gi)+t_div-t_gap));
                    
                    [sorted_cumLLR,sort_order] = sortrows(cumLLR_shuffled(pick_in_out & pick_num_accum,gi));
                    
                    sortedFirstShapeKM = oneShapeKM(pick_in_out & pick_num_accum,:);
                    sortedFirstShapeKM = sortedFirstShapeKM(sort_order,:);
                    
                    popID4KM = popID_shuffled(pick_in_out & pick_num_accum);
                    popID4KM = popID4KM(sort_order);

                    % create mean FR matrix for each cell
                    sortedFirstShapeKMmean = nan(size(sortedFirstShapeKM));
                    for k = 1:max(popID4KM)
                        pick = logical(popID4KM==k);
                        sortedFirstShapeKMmean(pick,:) = repmat(mean(double(sortedFirstShapeKM(pick,:)),1),sum(pick),1);
                    end
                    
                    pick = ~isnan(sorted_cumLLR);
                    sorted_cumLLR = sorted_cumLLR(pick);
                    sortedFirstShapeKM = sortedFirstShapeKM(pick,:);
                    num_per_group = floor(sum(pick)/group);
                    for gi = 1:group
                        pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                        woe_mean(i,gi) = mean(sorted_cumLLR(pick_group));
                        meanFirstShapeDeltaKM(i,gi,:) = mean(double(sortedFirstShapeKM(pick_group,:)),1)/FRgain;
                        sdFirstShapeDeltaKM(i,gi,:) = std(double(sortedFirstShapeKM(pick_group,:)),1)/FRgain;
                        seFirstShapeDeltaKM(i,gi,:) = sdFirstShapeDeltaKM(i,gi,:)./sqrt(num_per_group);
                        
                        % sdFirstShapeDeltaResKM(i,j,:) = std(double(sortedFirstShapeKM(pick_group,:)) - sortedFirstShapeKMmean(pick_group,:))/FRgain;
                        
                        epochGroupFR(i,gi,:) = mean(double(sortedFirstShapeKM(pick_group,t_pick)),2)/FRgain;
                        epochGroupFR_mean(i,gi) = mean(epochGroupFR(i,gi,:));
                        epochGroupFR_se(i,gi) = std(epochGroupFR(i,gi,:))/sqrt(num_per_group);
                        
                        % subplot(2,1,i)
                        % subplot(2,3,3*i-ei+1);
                        figure(fig);hold on;
                        plot_color = 1-(1-color_map(gi,:))/i^2; % Tout is plotted with a tint color.
                        fillTrace(t_axis(t_plot_pick),squeeze(meanFirstShapeDeltaKM(i,gi,t_div:t_div+250-t_gap)),squeeze(seFirstShapeDeltaKM(i,gi,t_div:t_div+250-t_gap)),plot_color);
                        % fillTrace(t_axis(t_plot_pick),squeeze(meanFirstShapeDeltaKM(i,j,t_div:t_div+250-t_gap)),squeeze(sdFirstShapeDeltaResKM(i,j,t_div:t_div+250-t_gap)),plot_color);
                        % plot(t_axis(t_plot_pick),squeeze(meanFirstShapeDeltaKM(i,j,t_div:t_div+250-t_gap)),'color',plot_color);
                        xlim([0 1000])
                    end
                    % Linear regression
                    [beta bint r rint stats] = regress(mean(sortedFirstShapeKM(:,t_pick),2)/FRgain,[ones(size(sorted_cumLLR)),sorted_cumLLR],0.05);
                    slope = beta(2)*10;
                    slopes_ci = (bint(2,2)-bint(2,1))/2*10;
                    offset = beta(1);
                    offset_ci = (bint(1,2)-bint(1,1))/2;
                    p_values = stats(3);
                    
                    FR_all_epoch{i} = [FR_all_epoch{i};mean(sortedFirstShapeKM(:,t_pick),2)/FRgain];
                    cumLLR_all_epoch{i} = [cumLLR_all_epoch{i};sorted_cumLLR];
                    epoch_label{i} = [epoch_label{i};gi*ones(size(sorted_cumLLR))];

                    figure(fig+1)
                    subplot(2,3,3*i-gi+1);hold on;
                    for gi = 1:group
                        h = ploterr(woe_mean(i,gi)/10,epochGroupFR_mean(i,gi),[],epochGroupFR_se(i,gi),1,'ko');
                        set(h(1),'MarkerSize',12,'color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:))
                    end
                    plot([woe_mean(i,1)/10,woe_mean(i,end)/10],beta(1)+slope*[woe_mean(i,1),woe_mean(i,end)]/10,'k-')
                    axis([-3 3 0 70]);
                    maxY = 70;
                    text(-2.5,maxY*(0.6*(i-1)+0.3),sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                    text(-2.5,maxY*(0.6*(i-1)+0.1),sprintf('%.1f\\pm%.1f',slope,slopes_ci),'FontSize',16,'FontWeight','bold');
                    text(0.5,maxY*(0.6*(i-1)+0.3),sprintf('Offset'),'FontSize',16,'FontWeight','bold');
                    text(0.5,maxY*(0.6*(i-1)+0.1),sprintf('%.1f\\pm%.1f',offset,offset_ci),'FontSize',16,'FontWeight','bold');
                    set(gca,'FontSize',16,'FontWeight','bold','Box','off','TickDir','out')
                   
                case 'LLR'
                    t_plot_pick = logical(t_axis>=t_axis(250*(3-gi)+t_div+10) & t_axis<=t_axis(250*(4-gi)+1+t_div-t_gap));
                    t_pick = logical(t_axis>=t_axis(t_div+10) & t_axis<=t_axis(250+1+t_div-t_gap));
                    
                    color_map = flipud(colormap(gray(10)));
                    color_map = color_map(2:end-1,:);
                    fig_h_size = 300;
                    switch i
                        case 1
                            woe2analyze = 1:8;
                            woe2regress = 5:8;
                            % woe2regress = 1:8;
                        case 2
                            % continue
                            woe2analyze = 1:8;
                            woe2regress = 1:4;
                            % woe2regress = 1:8;
                    end
                    color_ind = 1;
                    for gi = woe2analyze
                        pick_woe = logical(LLR_shuffled(:,gi)==woe(gi));
                        pick = logical(pick_in_out & pick_woe);
                        fr_picked = oneShapeKM(pick,:);                        
                        cumLLR_picked = subCumLLR_shuffled(pick,:);
                        
                        sac_fr_woe_mean(i,gi,:) = mean(double(fr_picked)/FRgain,1);
                        sac_fr_woe_sd(i,gi,:) = std(double(fr_picked)/FRgain,1);
                        sac_fr_woe_se(i,gi,:) = squeeze(sac_fr_woe_sd(i,gi,:))./sqrt(sum(pick));
                        
                        epochGroupDeltaFR = mean(double(fr_picked(:,t_pick))/FRgain,2) - mean(double(fr_picked(:,t_pre_pick))/FRgain,2);
                        epochGroupDeltaFR_mean(i,gi) = mean(epochGroupDeltaFR);
                        epochGroupDeltaFR_se(i,gi) = std(epochGroupDeltaFR)/sqrt(sum(pick));
                        
                        LLR_mean(i,gi,:) = nanmean(cumLLR_picked);
                        LLR_sd(i,gi,:) = nanstd(cumLLR_picked);
                        LLR_se(i,gi,:) = LLR_sd(i,gi,:)./sqrt(sum(pick));

                        lastStim2Sac_mean(i,gi) = mean(lastStim2Sac_shuffled(pick));
                        num_accum_group{gi,gi} = num_accum_shuffled(pick);
                        
                        figure(fig);hold on;
                        plot_color = 1-(1-color_map(color_ind,:))/i^2; % Tout is plotted with a tint color.
                        fillTrace(t_axis(t_plot_pick)',squeeze(sac_fr_woe_mean(i,gi,t_pick)),1*squeeze(sac_fr_woe_se(i,gi,t_pick)),plot_color);
                        % plot(t_axis(t_plot_pick)',squeeze(sac_fr_woe_mean(i,wi,t_pick)),'color',plot_color);
                        color_ind = color_ind+1;
                        xlim([0 1000])
                        % fillTrace(t_axis(t_plot_pick),squeeze(sac_fr_woe_mean(i,wi,t_div:t_div+250-t_gap)),0*squeeze(sac_fr_woe_se(i,wi,t_div:t_div+250-t_gap)),plot_color);
                    end
                    
                    pick = (pick_in_out & pick_sign);
                    % pick = pick_in_out;
                    epochDeltaFR = mean(double(oneShapeKM(pick,t_pick)),2)/FRgain... 
                                 - mean(double(oneShapeKM(pick,t_pre_pick)),2)/FRgain;
                    
                    % Linear regression
                    [beta bint r rint stats] = regress(epochDeltaFR,[ones(size(LLR_shuffled(pick,gi))),LLR_shuffled(pick,gi)],0.05);
                    
                    offset(gi,i) = beta(1);
                    offset_ci(gi,i) = (bint(1,2)-bint(1,1))/2;
                    slope(gi,i) = beta(2)*10;
                    slope_ci(gi,i) = (bint(2,2)-bint(2,1))/2*10;
                    p_values(gi,i) = stats(3);
                    
                    figure(fig+1)
                    % title(sprintf('%d-shape trial',ni),'FontSize',24,'FontWeight','bold');
                    subplot(2,3,3*i-gi+1);hold on;
                    color_ind = 1;
                    for gi = woe2analyze % 1:length(woe)
                        h = ploterr(woe(gi)/10,epochGroupDeltaFR_mean(i,gi),[],epochGroupDeltaFR_se(i,gi),1,'ko');
                        set(h(1),'MarkerSize',12,'color','k','MarkerFaceColor','k')
                        % set(h(1),'MarkerSize',12,'color',color_map(color_ind,:),'MarkerFaceColor',color_map(color_ind,:))
                        color_ind = color_ind+1;
                    end
                    plot([woe(woe2regress(1)),woe(woe2regress(end))]/10,offset(gi,i)+slope(gi,i)*[woe(woe2regress(1)),woe(woe2regress(end))]/10,'k-')
                    % axis([-1.1 1.1 -10 25]);
                    axis([-1.1 1.1 -25 25]);
                    maxY = 25;
                    text(-0.9,maxY*1,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                    text(-0.9,maxY*0.75,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
                    text(0.2,maxY*1,sprintf('%.1f\\pm%.1f',slope(gi,i),slope_ci(gi,i)),'FontSize',16,'FontWeight','bold');
                    text(0.2,maxY*0.75,sprintf('%.1f\\pm%.1f',offset(gi,i),offset_ci(gi,i)),'FontSize',16,'FontWeight','bold');
                    
                    set(gca,'FontSize',16,'FontWeight','bold','Box','off','TickDir','out')
            end
            
            figure(fig);hold on;
            switch i
                case 1
                    ylim([0 150])
                case 2
                    ylim([0 150])
            end
            switch gi
                case 1
                    xlabel('Time from N th shape onset','FontSize',24,'FontWeight','bold');
                otherwise
                    xlabel(sprintf('Time from N-%d th shape onset',gi-1),'FontSize',24,'FontWeight','bold'); 
            end
            % xlabel('Time from first shape onset','FontSize',24,'FontWeight','bold');
            ylabel('FR ± sem','FontSize',24,'FontWeight','bold');
            % title(sprintf('%d-shape trial',ni),'FontSize',24,'FontWeight','bold');
            set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
        end
    end
    
    figure(fig)
    set(gcf,'position',[100 600 fig_h_size 480],'Color','w')
    figure(fig+1)
    % supertitle(sprintf('%d-shape trial',ni));
    set(gcf,'position',[100 100 700 400],'Color','w')

    fig = fig+2;
% end
    
    if 0 % strcmp(sortBy,'cumLLR')
        % Testing whether FR convergence is more pronounced for Tin
        for i = 1:2
            [h,atab,ctab,stats] = aoctool(cumLLR_all_epoch{i}/10,FR_all_epoch{i},epoch_label{i},[],'cumLLR','FR','epoch');
        end
    end
    
    % plot subjective woe
    shape_order = [2,4,6,8,7,5,3,1];
    SWOE = shapeSWOE(shape_order);
    SWOE_se = shapeSWOE_se(shape_order);

    switch id
        case 1
            TinSWOE = flipud(-SWOE(1:4));
            TinSWOE_se = flipud(SWOE_se(1:4));
            ToutSWOE = flipud(-SWOE(5:8));
            ToutSWOE_se = flipud(SWOE_se(5:8));
        case 2
            TinSWOE = SWOE(5:8);
            TinSWOE_se = SWOE_se(5:8);
            ToutSWOE = SWOE(1:4);
            ToutSWOE_se = SWOE_se(1:4);
    end
    figure(fig+2)
    subplot(2,1,1);hold on;
    for gi = 1:4
        h=ploterr(woe(4+gi)/10,TinSWOE(gi),[],TinSWOE_se(gi),1,'o','abshhy',0);
        set(h(1),'Color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:),'MarkerSize',12)
        set(h(2),'Color','k')
    end
    xlim([0 1])
    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')

    subplot(2,1,2);hold on;
    for gi = 1:4
        h=ploterr(woe(gi)/10,ToutSWOE(gi),[],ToutSWOE_se(gi),1,'o','abshhy',0);
        set(h(1),'Color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:),'MarkerSize',12)
        set(h(2),'Color','k')
    end
    xlim([-1 0])
    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
    set(gcf,'position',[1400 100 400 800],'Color','w')
    
end

%% Causal saccadic FR conditioned on delta LLR
% ** I do not think this figure adds any extra support for the FR
% convergence at the end of decision...

if fig_switch(13)
    clearvars slopes slope_ci slopes_se p_values offsets offsets_ci
    clearvars meanFirstShapeDeltaKM sdFirstShapeDeltaKM seFirstShapeDeltaKM epochGroupFR epochGroupFR_mean epochGroupFR_se
                        
    num_fig = 3;
    FH = figure(fig);
    set(FH,'position',[550 100 400 800]);
    FH1=figure(fig+1);clf;hold on
    set(FH1,'position',[100 100 450 800]);

    num_trial = size(TM,1);
    group = 5;
    woe = [-9 -7 -5 -3 3 5 7 9];
    
    interval = floor((64-1)/(group-1));
    if ~exist('map','var')
        map = colormap;
    end
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    num_bin = (postsac_t-presac_t)/shift+1;
    plot_range = [200 -800; -850 -1050; -1200 -1400];
    % plot_range = [200 -450; -500 -700; -750 -950];
    % plot_range = [200 -400; -450 -650; -700 -900];
    % plot_range = [200 -350; -400 -600; -650 -850];
    % plot_range = [200 -300; -350 -550; -600 -800];
    % plot_range = [200 -250; -300 -500; -550 -750];
    % plot_range = [200 -200; -250 -450; -500 -700];
    
    sortBy = 'cumLLR';
    % Extract the peri-saccadic FR [-800ms 200ms]
    t_pre_sac = -800;
    t_post_sac = 400;
    t_axis = t_pre_sac:t_post_sac;
    t_gap = 20;
    sacKM = nan(size(KM,1),t_post_sac-t_pre_sac+1);
    for ti = 1:size(KM,1)
        RT_ind = round(RT(ti))+1;
        t_min_ind = RT_ind+t_pre_sac;
        t_max_ind = RT_ind+t_post_sac;
        if t_min_ind<1
            sacKM(ti,-t_min_ind+2:end) = KM(ti,1:t_max_ind);
        else
            sacKM(ti,:) = KM(ti,t_min_ind:t_max_ind);
        end
    end
    
    for gi = 1
        lastStim2Sac = TM(:,4);
        pick = logical(lastStim2Sac<lastAccumShape2Sac);
        lastStim2Sac(pick) = lastStim2Sac(pick)+250;
            
        % There is no point in shuffling, but still doing just in case...
        sacKM_shuffled = sacKM(rand_trial_order,:);
        T_in_out_shuffled = T_in_out(rand_trial_order);
        num_accum_shuffled = num_accum(rand_trial_order);
        popID_shuffled = popID(rand_trial_order');
        correct_shuffled = correct(rand_trial_order');
        LLR_shuffled = LLR{2}(rand_trial_order,:);
        cumLLR_shuffled = cumLLR{2}(rand_trial_order,:);
        subCumLLR_shuffled = subCumLLR{2}(rand_trial_order,:);
        lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
        
        for i = 1:2 % Tin or Tout
            clearvars epochGroupFR
            legend_str = [];
            switch i
                case 1
                    woe2analyze = 1:8;
                case 2
                    woe2analyze = 1:8;
            end
            pick_in_out = logical(T_in_out_shuffled==i); % changed!!
            pick_sign = logical((-1)^i*LLR_shuffled(:,gi)<0); % choose positive for Tin and negative for Tout
            t_pick = logical(t_axis>=-150 & t_axis<=-50);
            t_plot_pick = logical(t_axis>=-300 & t_axis<=200);
            switch sortBy
                case 'LLR'
                    color_map = (colormap(jet(8)));
                    % ******* plotting subjective WOE ******** 
    
                    shape_order = [2,4,6,8,7,5,3,1];
                    SWOE = shapeSWOE(shape_order);
                    SWOE_se = shapeSWOE_se(shape_order);

                    switch id
                        case 1
                            TinSWOE = flipud(-SWOE(1:4));
                            TinSWOE_se = flipud(SWOE_se(1:4));
                            ToutSWOE = flipud(-SWOE(5:8));
                            ToutSWOE_se = flipud(SWOE_se(5:8));
                        case 2
                            TinSWOE = SWOE(5:8);
                            TinSWOE_se = SWOE_se(5:8);
                            ToutSWOE = SWOE(1:4);
                            ToutSWOE_se = SWOE_se(1:4);
                    end
                    figure(fig)
                    subplot(2,1,1);hold on;
                    for gi = 1:4
                        h=ploterr(woe(4+gi)/10,TinSWOE(gi),[],TinSWOE_se(gi),1,'o','abshhy',0);
                        set(h(1),'Color',color_map(4+gi,:),'MarkerFaceColor',color_map(4+gi,:),'MarkerSize',12)
                        set(h(2),'Color','k')
                    end
                    xlim([0 1])
                    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
                    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')

                    subplot(2,1,2);hold on;
                    for gi = 1:4
                        h=ploterr(woe(gi)/10,ToutSWOE(gi),[],ToutSWOE_se(gi),1,'o','abshhy',0);
                        set(h(1),'Color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:),'MarkerSize',12)
                        set(h(2),'Color','k')
                    end
                    xlim([-1 0])
                    ylabel('Subjective WOE','FontSize',18,'FontWeight','bold');
                    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')

                    % ****************************************
                    
                    for gi = woe2analyze
                        pick_woe = logical(LLR_shuffled(:,gi)==woe(gi));
                        pick = logical(pick_in_out & pick_woe);
                        fr_picked = sacKM_shuffled(pick,:);
                        cumLLR_picked = subCumLLR_shuffled(pick,:);

                        n_sacKM = sum(isfinite(sacKM_shuffled(pick,:)));

                        sac_fr_woe_mean(i,gi,:) = nanmean(double(fr_picked)/FRgain);
                        sac_fr_woe_sd(i,gi,:) = nanstd(double(fr_picked)/FRgain);
                        sac_fr_woe_se(i,gi,:) = squeeze(sac_fr_woe_sd(i,gi,:))./sqrt(n_sacKM');
                        
                        % epochGroupDeltaFR = mean(double(fr_picked(:,t_pick))/FRgain,2) - mean(double(fr_picked(:,t_pre_pick))/FRgain,2);
                        % epochGroupDeltaFR_mean(i,wi) = mean(epochGroupDeltaFR);
                        % epochGroupDeltaFR_se(i,wi) = std(epochGroupDeltaFR)/sqrt(sum(pick));

                        LLR_mean(i,gi,:) = nanmean(cumLLR_picked);
                        LLR_sd(i,gi,:) = nanstd(cumLLR_picked);
                        LLR_se(i,gi,:) = LLR_sd(i,gi,:)./sqrt(sum(pick));

                        lastStim2Sac_mean(i,gi) = mean(lastStim2Sac_shuffled(pick));
                        
                        figure(fig+1)
                        subplot(2,1,i);hold on
                        fillTrace(t_pre_sac:t_post_sac,squeeze(sac_fr_woe_mean(i,gi,:)),0*squeeze(sac_fr_woe_se(i,gi,:)),color_map(gi,:));
                    end
                    
                case 'cumLLR'
                    fig_h_size = 800;
                    
                    % t_plot_pick = logical(t_axis>=t_axis(250*(3-ei)+t_div) & t_axis<=t_axis(250*(4-ei)+t_div-t_gap));
                    
                    [sorted_cumLLR,sort_order] = sortrows(cumLLR_shuffled(pick_in_out,gi));
                    sortedSacKM = sacKM_shuffled(pick_in_out,:);
                    sortedSacKM = sortedSacKM(sort_order,:);
                    pick = ~isnan(sorted_cumLLR);
                    sorted_cumLLR = sorted_cumLLR(pick);
                    sortedSacKM = sortedSacKM(pick,:);
                    num_per_group = floor(sum(pick)/group);
                    for gi = 1:group
                        pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                        woe_mean(i,gi) = mean(sorted_cumLLR(pick_group));
                        meanFirstShapeDeltaKM(i,gi,:) = mean(double(sortedSacKM(pick_group,:)),1)/FRgain;
                        sdFirstShapeDeltaKM(i,gi,:) = std(double(sortedSacKM(pick_group,:)),1)/FRgain;
                        seFirstShapeDeltaKM(i,gi,:) = sdFirstShapeDeltaKM(i,gi,:)./sqrt(num_per_group);
                        
                        epochGroupFR(i,gi,:) = mean(double(sortedSacKM(pick_group,t_pick)),2)/FRgain;
                        epochGroupFR_mean(i,gi) = mean(epochGroupFR(i,gi,:));
                        epochGroupFR_se(i,gi) = std(epochGroupFR(i,gi,:))/sqrt(num_per_group);
                        
                        % subplot(2,1,i)
                        % subplot(2,3,3*i-ei+1);
                        figure(fig+1);hold on;
                        plot_color = 1-(1-color_map(gi,:))/i^2; % Tout is plotted with a tint color.
                        fillTrace(t_axis(t_plot_pick),squeeze(meanFirstShapeDeltaKM(i,gi,t_plot_pick)),squeeze(seFirstShapeDeltaKM(i,gi,t_plot_pick)),plot_color);
                        xlim([-300 200])
                    end
                    % Linear regression
                    [beta bint r rint stats] = regress(mean(sortedSacKM(:,t_pick),2)/FRgain,[ones(size(sorted_cumLLR(:,gi))),sorted_cumLLR(:,gi)],0.05);
                    slope = beta(2)*10;
                    slope_ci = (bint(2,2)-bint(2,1))/2*10;
                    p_values = stats(3);

                    figure(fig+2)
                    subplot(2,3,3*i-gi+1);hold on;
                    for gi = 1:group
                        h = ploterr(woe_mean(i,gi)/10,epochGroupFR_mean(i,gi),[],epochGroupFR_se(i,gi),1,'ko');
                        set(h(1),'MarkerSize',12,'color',color_map(gi,:),'MarkerFaceColor',color_map(gi,:))
                    end
                    plot([woe_mean(i,1)/10,woe_mean(i,5)/10],beta(1)+slope*[woe_mean(i,1),woe_mean(i,5)]/10,'k-')
                    axis([-3 3 0 70]);
                    maxY = 70;
                    text(-2.5,maxY*(0.6*(i-1)+0.3),sprintf('Slope'),'FontSize',16,'FontWeight','bold');
                    text(-2.5,maxY*(0.6*(i-1)+0.1),sprintf('%.1f\\pm%.1f',slope,slope_ci),'FontSize',16,'FontWeight','bold');
                    set(gca,'FontSize',16,'FontWeight','bold','Box','off','TickDir','out')
            end
        end
    end
    figure(fig+1);
    set(gcf,'position',[100 100 500 600],'Color','w')
    
    fig = fig+num_fig;
    return
    for fi = 1:2
        for i = 1:2 % Tin or Tout
            color_ind = 2;
            switch i
                case 1
                    woe2analyze = 1:8;
                case 2
                    woe2analyze = 1:8;
            end
            for wi = woe2analyze
                switch fi
                    case 1
                        figure(fig+1)
                        subplot(2,1,i);hold on
                        fillTrace(t_pre_sac:t_post_sac,squeeze(sac_fr_woe_mean(i,wi,:)),0*squeeze(sac_fr_woe_se(i,wi,:)),color_map(wi,:));
                        % p=ploterr(t_pre_sac:t_post_sac,squeeze(sac_fr_woe_mean(i,wi,:)),[],1*squeeze(sac_fr_woe_se(i,wi,:)),1,'-','abshhy',0);
                        % pick = true(1,length(bin_center));
                        % set(p,'color',color_map(wi,:))
                        color_ind = color_ind+1;
                        hold on
                        plot([-lastStim2Sac_mean(i,wi) -lastStim2Sac_mean(i,wi)],[0 100],'--','color',color_map(wi,:))
                        
                        xlabel('Time from saccade','FontSize',18,'FontWeight','bold'); 
                        ylabel('FR ± SD','FontSize',18,'FontWeight','bold');
                        set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out')
                        xlim([-600 200])
                        switch i
                            case 1
                                ylim([20 70])
                            case 2
                                ylim([0 50])
                        end
                    case 2
                        figure(fig+3)
                        subplot(2,1,i);hold on
                        p = ploterr(1:2,flipud(squeeze(LLR_mean(i,wi,1:2))),[],flipud(squeeze(LLR_se(i,wi,1:2))),1,'-','abshhy',0);
                        set(p,'color',color_map(wi,:))
                        color_ind = color_ind+1;
                        xlim([0 3]);
                end
            end
        end
    end
    
    figure(fig+1);
    set(gcf,'position',[100 100 500 600],'Color','w')
    
    info.sac_fr_woe_mean = sac_fr_woe_mean;
    info.sac_fr_woe_sd = sac_fr_woe_sd;
    info.sac_fr_woe_se = sac_fr_woe_se;
    fig = fig+num_fig;
end

%% Causal saccadic FR vs cumulative WOE correlation (time stepped)
if fig_switch(14)
    clearvars slopes slope_ci slopes_se p_values offsets offsets_ci
    clearvars meanFirstShapeDeltaKM sdFirstShapeDeltaKM seFirstShapeDeltaKM epochGroupFR epochGroupFR_mean epochGroupFR_se
    
    num_fig = 1;
    figure(fig);hold on
    set(gcf,'position',[100 100 1000 500],'Color','w');

    group = 5;
    interval = floor((64-1)/(group-1));
    if ~exist('map','var')
        map = colormap;
    end
    color_map = map(1:interval:(1+interval*(group-1)),:);

    % Extract the peri-saccadic FR [-800ms 200ms]
    t_pre_sac = -800;
    t_post_sac = 400;
    t_axis = t_pre_sac:t_post_sac;

    bin_size = 20; % ms
    t_center_ind = 1:length(t_axis);
    % t_center_ind = (bin_size/2+1):(length(t_axis)-bin_size/2+1);
    t_center = t_axis(t_center_ind);
        
    if 1
        
        sacKM = nan(size(KM,1),t_post_sac-t_pre_sac+1);
        for ti = 1:size(KM,1)
            RT_ind = round(RT(ti))+1;
            t_min_ind = RT_ind+t_pre_sac;
            t_max_ind = RT_ind+t_post_sac;
            if t_min_ind<1
                sacKM(ti,-t_min_ind+2:end) = KM(ti,1:t_max_ind);
            else
                sacKM(ti,:) = KM(ti,t_min_ind:t_max_ind);
            end
        end

        lastStim2Sac = TM(:,4);

        % There is no point in shuffling, but still doing just in case...
        sacKM_shuffled = sacKM(rand_trial_order,:);
        T_in_out_shuffled = T_in_out(rand_trial_order);
        num_accum_shuffled = num_accum(rand_trial_order);
        popID_shuffled = popID(rand_trial_order');
        correct_shuffled = correct(rand_trial_order');
        rawCumLLR_shuffled = rawCumLLR{2}(rand_trial_order,:);
        LLR_shuffled = LLR{2}(rand_trial_order,:);
        cumLLR_shuffled = cumLLR{2}(rand_trial_order,:);
        subCumLLR_shuffled = subCumLLR{2}(rand_trial_order,:);
        lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
        RT_shuffled = RT(rand_trial_order);

        extra_delay = 100; %0:50:500;

        for di = 1:length(extra_delay)

            figure(fig);clf;hold on

            % create a matrix of the total logLR at each point in time.
            cumWoeKM = nan(size(KM,1),length(t_axis));
            sac_ind = -t_pre_sac + info.delay + extra_delay(di) + 1; % time when logLR at saccade is reflected in FR
            for ti = 1:size(KM,1)
                t_lastStim_ind = sac_ind-round(lastStim2Sac_shuffled(ti));
                if t_lastStim_ind<1
                    continue
                end
                cumWoeKM(ti,t_lastStim_ind:end) = rawCumLLR_shuffled(ti,1);
                for si = 1:5
                    if t_lastStim_ind-250*si>0
                        t_pick = (t_lastStim_ind-250*si):(t_lastStim_ind-250*(si-1)-1);
                        cumWoeKM(ti,t_pick) = rawCumLLR_shuffled(ti,si+1);
                    else
                        t_pick = 1:t_lastStim_ind-250*(si-1)-1;
                        cumWoeKM(ti,t_pick) = rawCumLLR_shuffled(ti,si+1);
                        break
                    end
                end
            end

            if 0 % plot example traces for cumWoe over time
                close all
                figure; hold on
                plot(cumWoeKM(1:20,:)')
                plot([1001,1001],[-30 30],'k--');
                plot([801,801],[-30 30],'r--');
                hold off
            end

            for i = [2,1] % Tin or Tout
                pick = logical(T_in_out_shuffled==i);
                % fr_sorted{i} = nan(sum(pick),length(t_center));
                % cumWoe_sorted{i} = nan(sum(pick),length(t_center));
                fr_mean_group{i} = nan(group,length(t_center));
                fr_se_group{i} = nan(group,length(t_center));

                for ti = 1:length(t_center);
                    t_pick = ti;
                    % t_pick = (t_center_ind(ti)-bin_size/2):(t_center_ind(ti)+bin_size/2-1);

                    cumWoe_pick = mean(cumWoeKM(pick,t_pick),2);
                    fr_pick = mean(sacKM_shuffled(pick,t_pick),2);

                    pick_finite = isfinite(cumWoe_pick);
                    cumWoe_finite = cumWoe_pick(pick_finite);
                    fr_finite = fr_pick(pick_finite);

                    [cumWoe_sorted,sort_order] = sortrows(cumWoe_finite);
                    fr_sorted = fr_finite(sort_order);

                    num_per_group = floor(sum(pick_finite)/group);
                    for gi = 1:group
                        pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                        fr_group{i,gi,ti} = fr_sorted(pick_group);
                        fr_mean_group{i}(gi,ti) = nanmean(fr_sorted(pick_group)/FRgain);
                        fr_se_group{i}(gi,ti) = nanstd(fr_sorted(pick_group)/FRgain)/sqrt(num_per_group);
                    end

                    % [cumWoe_sorted{i}(:,ti),sort_order] = sortrows(cumWoe_pick);
                    % fr_sorted{i}(:,ti) = fr_pick(sort_order);
                end
                for gi = 1:group
                    plot_color = 1-(1-color_map(gi,:))/i^2; % Tout is plotted with a tint color.
                    fillTrace(t_center,fr_mean_group{i}(gi,:),fr_se_group{i}(gi,:),plot_color);
                end
            end

            hold on
            plot([-200 -200],[0 100],'k--');
            plot([0 0],[0 100],'k--');
            ylim([0 70])
            text(200,70*0.9,[num2str(info.delay+extra_delay(di)),'ms lag'],'FontSize',24,'FontWeight','bold')

            % save the figure in eps format
            % figpath = '/Users/shin/Documents/RTshape/RTshapeFigs2/';
            % print2eps([figpath,'Joey_sacFR_lag_',num2str(info.delay+extra_delay(di)),'ms.eps'],gcf)

            % save([save_dir,'EliSacFR'],'fr_group');
        end
        
    else    
        % combine FR from two monkeys
        load EliSacFR
        fr_group_1 = fr_group;
        load JoeySacFR
        fr_group_2 = fr_group;

        clearvars fr_group
        group = 5;
        for i = [2,1]
            for gi = 1:group
                for ti = 1:length(t_center);
                    if 0
                        % pool trials from two monkeys
                        fr_group_comb = [fr_group_1{i,gi,ti};fr_group_2{i,gi,ti}];
                        fr_mean_group{i}(gi,ti) = nanmean(fr_group_comb/FRgain);
                        fr_se_group{i}(gi,ti) = nanse(fr_group_comb/FRgain);
                    else
                        % simply average the mean FR from two monkeys
                        fr_mean_group{i}(gi,ti) = mean([nanmean(fr_group_1{i,gi,ti}),nanmean(fr_group_2{i,gi,ti})]/FRgain);
                        fr_se_group{i}(gi,ti) = 0;
                    end
                end
                plot_color = 1-(1-color_map(gi,:))/i^2; % Tout is plotted with a tint color.
                fillTrace(t_center,fr_mean_group{i}(gi,:),fr_se_group{i}(gi,:),plot_color);hold on;
            end
        end
        plot([-200 -200],[0 100],'k--');
        plot([0 0],[0 100],'k--');
        ylim([0 70])
    
    end
    fig = fig+num_fig;

end

%% Causal saccadic FR sorted by RT
if fig_switch(15)
    clearvars slopes slope_ci slopes_se p_values offsets offsets_ci
    clearvars meanFirstShapeDeltaKM sdFirstShapeDeltaKM seFirstShapeDeltaKM epochGroupFR epochGroupFR_mean epochGroupFR_se
    
    num_fig = 2;
    figure(fig);hold on
    set(gcf,'position',[100 100 1000 500],'Color','w');
    
    % Extract the peri-saccadic FR [-800ms 200ms]
    t_pre_sac = -800;
    t_post_sac = 400;
    t_axis = t_pre_sac:t_post_sac;

    bin_size = 1; % ms
    t_center_ind = 1:length(t_axis);
    t_center = t_axis(t_center_ind);
    
    sacKM = nan(size(KM,1),t_post_sac-t_pre_sac+1);
    for ti = 1:size(KM,1)
        RT_ind = round(RT(ti))+1;
        t_min_ind = RT_ind+t_pre_sac;
        t_max_ind = RT_ind+t_post_sac;
        if t_min_ind<1
            sacKM(ti,-t_min_ind+2:end) = KM(ti,1:t_max_ind);
        else
            sacKM(ti,:) = KM(ti,t_min_ind:t_max_ind);
        end
    end

    lastStim2Sac = TM(:,4);

    % There is no point in shuffling, but still doing just in case...
    sacKM_shuffled = sacKM(rand_trial_order,:);
    T_in_out_shuffled = T_in_out(rand_trial_order);
    num_accum_shuffled = num_accum(rand_trial_order);
    popID_shuffled = popID(rand_trial_order');
    correct_shuffled = correct(rand_trial_order');
    rawCumLLR_shuffled = rawCumLLR{2}(rand_trial_order,:);
    LLR_shuffled = LLR{2}(rand_trial_order,:);
    cumLLR_shuffled = cumLLR{2}(rand_trial_order,:);
    subCumLLR_shuffled = subCumLLR{2}(rand_trial_order,:);
    lastStim2Sac_shuffled = lastStim2Sac(rand_trial_order);
    RT_shuffled = RT(rand_trial_order);

    figure(fig);clf;hold on
    
    switch id
        case 1
            ei_set = 8:-1:3;
        case 2
            ei_set = 6:-1:1;
    end
    
    group = length(ei_set);
    interval = floor((64-1)/(group-1));
    if ~exist('map','var')
        map = colormap;
    end
    color_map = flipud(map(1:interval:(1+interval*(group-1)),:));
    
    extra_delay = 0;
    cumWoeKM = nan(size(KM,1),length(t_axis));
    sac_ind = -t_pre_sac + info.delay + extra_delay + 1; % time when logLR at saccade is reflected in FR
    for ti = 1:size(KM,1)
        t_lastStim_ind = sac_ind-round(lastStim2Sac_shuffled(ti));
        if t_lastStim_ind<1
            continue
        end
        cumWoeKM(ti,t_lastStim_ind:end) = rawCumLLR_shuffled(ti,1);
        for si = 1:5
            if t_lastStim_ind-250*si>0
                t_pick = (t_lastStim_ind-250*si):(t_lastStim_ind-250*(si-1)-1);
                cumWoeKM(ti,t_pick) = rawCumLLR_shuffled(ti,si+1);
            else
                t_pick = 1:t_lastStim_ind-250*(si-1)-1;
                cumWoeKM(ti,t_pick) = rawCumLLR_shuffled(ti,si+1);
                break
            end
        end
    end
    sacKM_shuffled(isnan(cumWoeKM)) = nan;
    
    for i = [2,1] % Tin or Tout
        for gi = 1:length(ei_set)
            pick = logical(T_in_out_shuffled==i & num_accum_shuffled==ei_set(gi));
            fr_mean_group{i}(gi,:) = nanmean(sacKM_shuffled(pick,:)/FRgain);
            fr_se_group{i}(gi,:) = 1*nanse(sacKM_shuffled(pick,:)/FRgain);
            
            t_plot_ind = logical(sum(isfinite(sacKM_shuffled(pick,:)))>(sum(pick)/3));
            
            % Do not show PSTH before FR starts to reflect the 1st shape's logLR
            % if (lastAccumShape2Sac + 250*ei_set(gi) - info.delay) <= -t_pre_sac
            %     t_ind = sac_ind-(lastAccumShape2Sac + 250*ei_set(gi) - info.delay); % Index of time when 1st shape starts to affect FR.
            %     fr_mean_group{i}(gi,1:t_ind) = nan;
            %     fr_se_group{i}(gi,1:t_ind) = nan;
            % end
            
            plot_color = 1-(1-color_map(gi,:))/i^2; % Tout is plotted with a tint color.
            fillTrace(t_center(t_plot_ind),fr_mean_group{i}(gi,t_plot_ind),fr_se_group{i}(gi,t_plot_ind),plot_color);
            hold on
        end
    end
   
    hold on
    plot([-200 -200],[0 100],'k--');
    plot([0 0],[0 100],'k--');
    ylim([0 80])
%%
    % create legend
    figure(fig+1);hold on
    for gi = 1:length(ei_set)
        plot_color = 1-(1-color_map(gi,:)); % Tout is plotted with a tint color.
        plot([0 1],[0 1],'-','color',plot_color,'LineWidth',3);
    end
    % legend('N*=3','N*=4','N*=5','N*=6','N*=7','N*=8','Location','NorthEastOutside');
    legend('N*=1','N*=2','N*=3','N*=4','N*=5','N*=6','Location','NorthEastOutside');
    
    % save the figure in eps format
    % figpath = '/Users/shin/Documents/RTshape/RTshapeFigs2/';
    % print2eps([figpath,'Joey_sacFR_lag_',num2str(info.delay+extra_delay(di)),'ms.eps'],gcf)

    % save([save_dir,'EliSacFR'],'fr_group');
        
    fig = fig+num_fig;

end


%% FR vs totalWOE partial correlation
if fig_switch(16)
    num_fig = 1;
    dw = 100; % ms
    dw_shift = 50; % ms
    t_axis = 0:dw_shift:dw_shift*100;
    figure(info.fig);clf;
    plot_count = 1;
    for num2pick = 3:8  
        for i = 1:num2pick+1;
            pick = find(num_accum==num2pick);
            fr_pick = FM(pick,:);
            cumWoe_pick = rawCumLLR{1}(pick,1:num2pick+1);
            k = setdiff(1:num2pick+1,i);
            control = cumWoe_pick(:,k);
            R = partialcorr(cumWoe_pick,fr_pick(:,1:50),control,'rows','pairwise');
            % R = partialcorr([cumWoe_pick,fr_pick(:,1:50)],'rows','pairwise');
            PM(:) = R(i,:);
        end
        subplot(2,3,plot_count);hold on
        imagesc(t_axis(1:50),1:num2pick+1,PM(:,1:50),[-0.6 0.6]);
        for i = 1:num2pick+1
            plot((250*i):dw_shift:(250*(i+1)-dw_shift),i,'ko','MarkerFaceColor','k');
        end
        sac_line = 250*(num2pick)+150;
        plot([sac_line,sac_line],[0 10],'k--');
        plot([sac_line+250,sac_line+250],[0 10],'k--');
        axis([-50 2500 0 10])
        plot_count = plot_count+1;
        info.PM = PM;
    end
    fig = fig+num_fig;
end

%% VarCE analysis
if fig_switch(17)
    num_fig = 4;
    fit_info.file_name = 'fitPhi';
    fit_info.file_path = save_dir;
    cell_ID = unique(popID);
    MinTrials = 10; % replace NaNs with counts only if there are at least this many trials
    
    % Calculate phi from FR around the target onset
    spTarg = TM(:,51)*0.2;
    
    switch id
        case 1
            minNstar = 10;
            NstarV = 4:8;
            BV = 14;
            this_target = 1; % left target (Tin)
        case 2
            minNstar = 4;
            NstarV = 1:5;
            BV = 8;
            this_target = 2; % right target (Tin)
    end
 
    for i_cell = 1:length(cell_ID)
        pick = logical(popID==cell_ID(i_cell));
        if sum(pick) >= MinTrials
            q = spTarg(pick,:);
            phi(i_cell) = nanvar(q)/nanmean(q);
            q_n(i_cell) = sum(isfinite(q)); % count
        end
    end
    % compute phi by weighted average
    phiTarg = q_n/sum(q_n)*phi';
    
    %%%%%%%%%% Forward analysis to fit Phi %%%%%%%%%
    spMtrx = frMtrx{1}*250/1e3;
    
    if 1
        for k = 1:size(spMtrx,1)
            spMtrx(k,num_accum(k):end) = nan;
        end
    end
    
    % Preallocation of memory
    mRaw = nan(size(spMtrx));
    mMean = nan(size(spMtrx));
    N = nan(size(spMtrx));
    
    for i_cell = 1:length(cell_ID)
        pick = logical(popID==cell_ID(i_cell) & num_accum>=1 & rew_targ==this_target);
        if sum(pick) >= MinTrials
            Q = spMtrx(pick,:);
            Qmean = repmat(nanmean(Q),size(Q,1),1); % mean -> nanmean? SK
            Qn = repmat(sum(isfinite(Q)),size(Q,1),1); % count
            % place into large area
            mRaw(pick,:) = Q;
            mMean(pick,:) = Qmean;
            N(pick,:) = Qn;
        end
    end
    
    pick = logical(rew_targ==this_target);
    xRaw = mRaw(pick,1:minNstar);
    xMean = mMean(pick,1:minNstar);
    xN = N(pick,1:minNstar);
    phi = nanvar(mRaw(pick,end))/nanmean(mRaw(pick,end));
    % C = VarCorCEshapeMixture(BV,0);
    % C = CorCEpartialSums(10);
    C = VarCorCEshape([14,3]);
    predCorCE = C(1:minNstar,1:minNstar);
    phiGrid = (.3:.03:.7)'; % to find a good starting guess
    if 1
        f = @(phi) bCorVarCEerr(phi,xRaw,xMean,xN,predCorCE);

        for i = 1:length(phiGrid)  % sample phi vals in a grid
            errG(i) = f(phiGrid(i));
        end
        [errGmin J] = min(errG);
        phiGbest = phiGrid(J);
        [phiFit,errFit] = fminsearch(f,phiGbest); % use the best to start fminsearch
    elseif 0
        C = 0; % Resetting the iteration counter
        Err = 0; % Resetting the error
        Prmts = zeros(1,3);
        
        save([fit_info.file_path,fit_info.file_name,'.mat'],'C', 'Prmts', 'Err');

        thetaGuess = [1,10,2];
        f = @(theta) bCorVarCEerrShape([theta(1),12,3],xRaw,xMean,xN,predCorCE,fit_info);

        [thetaFit,errFit] = fminsearch(f,thetaGuess);
        save([fit_info.file_path,fit_info.file_name,'.mat'],'thetaFit', 'errFit','-append');
        phi = thetaFit(1);
        % if errFit > errGmin
        %     warning('grid outperformed search')
        %     phiFit = phiGbest;
        %     errFit = errGmin;
        % end
    elseif 0
        load([fit_info.file_path,fit_info.file_name,'.mat'])
        thetaFit = Prmts(end,:);
        % phiFit = thetaFit(1);
        phiFit = 1;
    end
    phiFit = 0.55;
    [VarCEfit,CovCEfit,CorCEfit] = VarCovCorCE(xRaw,xMean,phiFit);
    fig = fig+2;
    h1=figure(fig);
    set(h1,'position',[0 50 400 800])
    subplot(2,1,1)
    imagesc(CorCEfit,[0 1]);
    colorbar;
    
    font_size = 16;
    set(gca,'XTick',1:minNstar,'XTickMode','manual','XTickLabel',1:minNstar);
    set(gca,'YTick',1:minNstar,'YTickMode','manual','YTickLabel',1:minNstar)
    title(sprintf('CorCE   phi = %.2f',phiFit),'FontSize',font_size,'FontWeight','bold')
    xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    ylabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')
    
    subplot(2,1,2)
    plot(1:minNstar,VarCEfit,'-k')
    xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    ylabel('VarCE','FontSize',font_size,'FontWeight','bold')
    set(gca,'XTick',1:minNstar,'XTickMode','manual','XTickLabel',1:minNstar);
    set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')
    shg
    
    %%%%%%%%% Backward analysis %%%%%%%%%%%%
    t_align = 0;% -200; CHANGED!!
    t_center = t_align-250*[0:10];
    for i = 1:length(t_center)
        t_ind(i) = find(info.sac_bin_center==t_center(i));
    end
    spMtrx = sacMtrx(:,t_ind)*info.binsize/1e3;
    
    % Preallocation of memory
    mRaw = nan(size(spMtrx));
    mMean = nan(size(spMtrx));
    N = nan(size(spMtrx));
    
    cell_ID = unique(popID);
    MinTrials = 10; % replace NaNs with counts only if there are at least this many trials
    phiGrid = (.3:.03:.7)'; % to find a good starting guess
    for i_cell = 1:length(cell_ID)
        for i_IO = 1:2
            for i_N = 1:8
                % pick = logical(popID==cell_ID(i_cell) & rew_targ==this_target);
                pick = logical(popID==cell_ID(i_cell) & rew_targ==this_target & T_in_out==i_IO & num_accum==i_N);
                if sum(pick) >= MinTrials
                    Q = spMtrx(pick,:);
                    Qmean = repmat(nanmean(Q),size(Q,1),1); % mean -> nanmean? SK
                    Qn = repmat(sum(isfinite(Q)),size(Q,1),1); % count
                    % place into large area
                    mRaw(pick,:) = Q;
                    mMean(pick,:) = Qmean;
                    N(pick,:) = Qn;
                end
            end
        end
    end
    
    VarCE_Mtrx = nan(length(NstarV),NstarV(end));
    CorCE_Mtrx = nan(NstarV(end),NstarV(end),length(NstarV));
    
    ind = 1;
    fig = fig+1;
    h1=figure(fig);
    set(h1,'position',[0*ind 200 1600 300])
    % Plot the diffusion process by including Tin & Tout
    subplot(1,length(NstarV)+1,ind)
    imagesc(CorCEfit,[0 1]);
    CorCE_Mtrx(1:minNstar,1:minNstar,1) = CorCEfit;
    axis square
    ax = gca;

    %set(gca,'XTick',1:minNstar,'XTickMode','manual','XTickLabel',1:minNstar);
    %set(gca,'YTick',1:minNstar,'YTickMode','manual','YTickLabel',1:minNstar)
    title(sprintf('CorCE   phi = %.2f',phiFit),'FontSize',font_size,'FontWeight','bold')
    xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    ylabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')
    
    sTIO = {'Tin','Tout','Tin & Tout'};
    % Plot the bounded accumulation by plotting Tin or Tout only
    for i_IO = 1
        for i_N = NstarV
            for i_cell = 1:length(cell_ID)
                if i_IO<3
                    pick = logical(rew_targ==this_target & T_in_out==i_IO & num_accum == i_N);
                else
                    pick = logical(rew_targ==this_target & num_accum==i_N);
                end
                
            end
            % pick = logical(num_accum == i_N);
            % pick = logical(T_in_out==i_IO & num_accum == i_N);
            xRaw = fliplr(mRaw(pick,1:i_N));
            xMean = fliplr(mMean(pick,1:i_N));
            xN = fliplr(N(pick,1:i_N));
            [VarCEfit,CovCEfit,CorCEfit] = VarCovCorCE(xRaw,xMean,phiFit);
            nanmean(xRaw)
            nanmean(xMean)
            VarCE_Mtrx(ind,NstarV(end)-i_N+1:NstarV(end)) = VarCEfit;
            CorCE_Mtrx(1:i_N,1:i_N,ind+1) = CorCEfit;
            ind = ind+1;
            subplot(1,length(NstarV)+1,ind);
            if isreal(CorCEfit)
                imagesc(CorCEfit,[0 1]);
            end
            axis square
            
            % title(sprintf('CorCE   phi = %.2f',phiFit),'FontSize',font_size,'FontWeight','bold')
            title(sprintf([sTIO{i_IO},'  N* = %d'],i_N),'FontSize',font_size,'FontWeight','bold')
            xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
            % ylabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
            set(gca,'XTick',1:i_N,'XTickMode','manual','XTickLabel',1:i_N);
            set(gca,'YTick',1:i_N,'YTickMode','manual','YTickLabel',1:i_N)
            set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')
            
            % subplot(2,1,2)
            % plot(1:i_N,VarCEfit,'-k')
            % xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
            % ylabel('VarCE','FontSize',font_size,'FontWeight','bold')
            % set(gca,'XTick',1:i_N,'XTickMode','manual','XTickLabel',1:i_N);
            % set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')
        end
    end

    fig = fig+1;
    figure(fig)
    plot(1:NstarV(end),VarCE_Mtrx','LineWidth',3);
    xtickStr{1} = 'N*';
    for i = 1:NstarV(end)-1
        xtickStr{i+1} = ['N*-',num2str(i)];
    end
    xlabel('Shape epoch #','FontSize',font_size,'FontWeight','bold')
    ylabel('VarCE','FontSize',font_size,'FontWeight','bold')
    set(gca,'XTick',1:NstarV(end),'XTickMode','manual','XTickLabel',fliplr(xtickStr));
    set(gca,'FontSize',font_size,'FontWeight','bold','Box','OFF','TickDir','out')

    fig = fig+1;
end

%% Forward mixture vs non-mixture (random walk) test
% Excluding the last shape to remove the saccadic contamination
if fig_switch(18)
    num_fig = 2;
    FH1 = figure(fig); clf;
    set(FH1,'position',[100 100 1000 400]);
    FH2 = figure(fig+1); clf; hold on;
    set(FH2,'position',[150 150 1000 400]);
    FH3 = figure(fig+2); clf; hold on;
    set(FH3,'position',[200 200 1000 400]);
    FH4 = figure(fig+3); clf; hold on;
    set(FH4,'position',[50 50 1600 1200]);
    
    
    uniqueNum = unique(num_accum);
    group = 5; % Set to 0 to plot for all unique logLR

    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    uniqueID = unique(info.popID);
    num2analyze = 6;
    fr_bin = -10:5:70;
    
    figure(fig)
    for ci = 1% :length(uniqueID)
        clearvars woe_pick woe_mean woe_var fr_pick fr_hist fr_woe_mean fr_woe_var
        % pick_trial = true(length(popID),1);
        % pick_trial = logical(num_accum==num2analyze);
        % pick_trial = logical(num_accum==num2analyze);
        pick_trial = logical(rew_targ==1 & num_accum==num2analyze);
        % pick_trial = logical(rew_targ==1 & correct & popID==uniqueID(ci) & num_accum==num2analyze);
        for gi = 1:num2analyze %uniqueNum(end)
            frMtrxTemp = [frMtrx{1}(pick_trial,gi),cumLLR{1}(pick_trial,gi)];
            pick = logical(~isnan(frMtrxTemp(:,end)) & frMtrxTemp(:,end)>-9000 & frMtrxTemp(:,end)<9000);
            % pick = logical(rew_targ(pick_trial)==1 & ~isnan(frMtrxTemp(:,end)) & frMtrxTemp(:,end)>-9000 & frMtrxTemp(:,end)<9000);
            sortedFrMtrx = sortrows(frMtrxTemp(pick,:),2);
            num_per_group = floor(sum(pick)/group);
            for gi = 1:group
                pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                woe_pick{gi,gi} = sortedFrMtrx(pick_group,2);
                woe_mean(gi,gi) = nanmean(woe_pick{gi,gi})./10;
                woe_var(gi,gi) = nanvar(woe_pick{gi,gi})./10;
                fr_pick{gi,gi} = sortedFrMtrx(pick_group,1);
                fr_hist(gi,:,gi) = histc(fr_pick{gi,gi},fr_bin);
                fr_woe_mean(gi,gi) = nanmean(fr_pick{gi,gi});
                fr_woe_var(gi,gi) = nanvar(fr_pick{gi,gi});
            end
            norm_fr_mean(gi,:,ci) = fr_woe_mean(gi,:)/mean(fr_woe_mean(gi,:));
            norm_fr_var(gi,:,ci) = fr_woe_var(gi,:)/mean(fr_woe_var(gi,:));
            subplot(2,4,gi);
            for k = 1:group
                plot(norm_fr_mean(gi,k,ci),norm_fr_var(gi,k,ci),'ko','MarkerSize',8,'MarkerFaceColor',color_map(k,:));
                hold on
                if k==1
                    xlabel('Normalized mean','FontSize',24,'FontWeight','bold'); 
                    ylabel('Normalized var','FontSize',24,'FontWeight','bold');
                end
            end
            if gi==1
                x_range = get(gca,'Xlim');
                y_range = get(gca,'Ylim');
            end
            % xlim([x_range(1)-10,x_range(2)+30])
            % ylim([0,y_range(2)+100])
            subplot(2,4,gi);hold off
        end
    end
    woe_pop_mean = nanmean(woe_mean,3);
    fr_pop_mean = nanmean(norm_fr_mean,3);
    fr_pop_var = nanmean(norm_fr_var,3);
    
    figure(fig+1)
    for gi = 1:num2analyze
        for i = 1:group
            subplot(2,4,gi); hold on
            plot(fr_pop_mean(gi,i),fr_pop_var(gi,i),'ko','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
            % xlim([0 2])
        end
    end
    
    figure(fig+2)
    for gi = 1:num2analyze
        for i = 1:group
            subplot(2,4,gi); hold on
            plot(woe_pop_mean(gi,i),fr_pop_mean(gi,i),'ko','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
            % xlim([-2 2]);
            % ylim([0 2])
        end
    end

    figure(fig+3)
    for gi = 1:num2analyze
        subplot(2,4,gi);
        h=bar3(fr_bin,fr_hist(:,:,gi)',1);
        if 1
            shading interp
            for i = 1:length(h)
                zdata = get(h(i),'ZData');
                set(h(i),'CData',zdata)
                % Add back edge color removed by interpolating shading
                set(h,'EdgeColor','k') 
            end
        end
        axis xy
        axis square
        view(30,80)
    end
    AZ = -40;
    EL = 45;
    if 0
        for ai = 1:36;
            AZ = AZ+10;
            for ei = 1:num2analyze
                figure(fig+3)
                subplot(2,4,ei)
                view(AZ,EL)
                drawnow
            end
            % pause(0.1)
        end
    end

    fig = fig+num_fig;
    
end
%% Backward mixture vs non-mixture (random walk) test
% Excluding the last shape to remove the saccadic contamination
if fig_switch(19)
    num_fig = 2;
    FH1 = figure(fig); clf;
    set(FH1,'position',[100 100 1000 400]);
    FH2 = figure(fig+1); clf; hold on;
    set(FH2,'position',[150 150 1000 400]);
    FH3 = figure(fig+2); clf; hold on;
    set(FH3,'position',[200 200 1000 400]);
    FH4 = figure(fig+3); clf; hold on;
    set(FH4,'position',[250 250 1000 400]);
    
    uniqueNum = unique(num_accum);
    group = 5; % Set to 0 to plot for all unique logLR

    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    uniqueID = unique(info.popID);
    num2analyze = 4;
    
    figure(fig)
    for ci = 1:length(uniqueID)
        % pick = true(length(popID),1);
        pick = logical(popID==uniqueID(ci));
        for gi = 1:num2analyze %uniqueNum(end)
            frMtrxTemp = [frMtrx{2}(pick,gi),cumLLR{2}(pick,gi)];
            pick = logical(~isnan(frMtrxTemp(:,end)) & frMtrxTemp(:,end)>-9000 & frMtrxTemp(:,end)<9000);
            sortedFrMtrx = sortrows(frMtrxTemp(pick,:),2);
            num_per_group = floor(sum(pick)/group);
            for gi = 1:group
                pick_group = (num_per_group*(gi-1)+1):(num_per_group*gi);
                woe_pick{gi,gi} = sortedFrMtrx(pick_group,2);
                woe_mean(gi,gi) = nanmean(woe_pick{gi,gi})./10;
                woe_var(gi,gi) = nanvar(woe_pick{gi,gi})./10;
                fr_pick{gi,gi} = sortedFrMtrx(pick_group,1);
                fr_hist(gi,:,gi) = histc(fr_pick{gi,gi},-10:2:70);
                fr_woe_mean(gi,gi) = nanmean(fr_pick{gi,gi});
                fr_woe_var(gi,gi) = nanvar(fr_pick{gi,gi});
            end
            norm_fr_mean(gi,:,ci) = fr_woe_mean(gi,:);
            norm_fr_var(gi,:,ci) = fr_woe_var(gi,:);
            norm_fr_mean(gi,:,ci) = fr_woe_mean(gi,:)./mean(fr_woe_mean(gi,:));
            norm_fr_var(gi,:,ci) = fr_woe_var(gi,:)./mean(fr_woe_var(gi,:));
            subplot(2,4,gi);hold on
            for k = 1:group
                plot(norm_fr_mean(gi,k,ci),norm_fr_var(gi,k,ci),'ko','MarkerSize',8,'MarkerFaceColor',color_map(k,:));
                if k==1
                    plot(norm_fr_mean(gi,:,ci),norm_fr_var(gi,:,ci),'k-','MarkerSize',8,'MarkerFaceColor',color_map(k,:));
                    xlabel('Normalized mean','FontSize',18,'FontWeight','bold'); 
                    ylabel('Normalized var','FontSize',18,'FontWeight','bold');
                end
            end
            if gi==1
                x_range = get(gca,'Xlim');
                y_range = get(gca,'Ylim');
            end
            xlim(x_range)
            ylim([0,y_range(2)])
        end
    end
    
    woe_pop_mean = nanmean(woe_mean,3);
    woe_pop_var = nanmean(woe_var,3);
    fr_pop_mean = nanmean(norm_fr_mean,3);
    fr_pop_var = nanmean(norm_fr_var,3);
    figure(fig+1)
    for gi = 1:num2analyze
        for i = 1:group
            subplot(2,4,gi); hold on
            var_mixed(gi,i) = mixtureVar(fr_pop_mean(gi,i),fr_pop_mean(gi,1),fr_pop_mean(gi,5),fr_pop_var(gi,1),fr_pop_var(gi,5));
            plot(fr_pop_mean(gi,i),fr_pop_var(gi,i),'ko','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
            plot(fr_pop_mean(gi,i),var_mixed(gi,i),'k^','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
            % xlim([0 2])
            % ylim([0 3])
        end
        if gi==1
            x_range = get(gca,'Xlim');
            y_range = get(gca,'Ylim');
        end
        % xlim(x_range)
        % ylim([0,y_range(2)])
        xlabel('Normalized FR mean','FontSize',14,'FontWeight','bold'); 
        ylabel('Normalized FR var','FontSize',14,'FontWeight','bold');
    end
    
    figure(fig+2)
    for gi = 1:num2analyze
        for i = 1:group
            subplot(2,4,gi); hold on
            plot(woe_pop_mean(gi,i),fr_pop_mean(gi,i),'ko','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
            xlim([-3 3]);
            ylim([0 2]);
        end
        xlabel('total logLR','FontSize',14,'FontWeight','bold'); 
        ylabel('Normalized FR mean','FontSize',14,'FontWeight','bold');
    end
    
    figure(fig+3)
    for gi = 1:num2analyze
        subplot(2,4,gi); hold on
        imagesc(1:5,-10:2:70,fr_hist(:,:,gi)')
        [errH0,errH1] = bimodalFtest(fr_pick{gi,3},fr_pick{gi,1},fr_pick{gi,5});
        Err(gi,:) = [errH0,errH1];
    end
    
%     figure(fig+3)
%     for ei = 1:num2analyze
%         for i = 1:group
%             subplot(2,4,ei); hold on
%             plot(fr_pop_mean(ei,i),woe_pop_var(ei,i),'ko','MarkerSize',8,'MarkerFaceColor',color_map(i,:))
%             xlim([0 2])
%         end
%         xlabel('total logLR','FontSize',14,'FontWeight','bold'); 
%         ylabel('logLR var','FontSize',14,'FontWeight','bold');
%     end

    fig = fig+num_fig;
    
end

%% Epoch-wise frequency histogram as a function of WOE & FR
if fig_switch(20)
    FH1 = figure(fig);
    set(FH1,'position',[100 100 1500 300])
    FH2 = figure(fig+1);
    set(FH2,'position',[100 200 1500 300])
    TinRew = logical(TM(:,1)==TM(:,7));
    ToutRew = logical(TM(:,1)~=TM(:,7));
    AUC_woe = nan(1,10);
    AUC_fr  = nan(1,10);
    frMtrx{1} = round(frMtrx{1}); % changed!
    cor_woe_all = [];
    cor_fr_all = [];
    res_woe_all = [];
    res_fr_all = [];
    for gi = 1:10
        woe_conv = -9*gi:2:9*gi;
        fr_conv = 0:1:100;
        zero_temp = zeros(length(fr_conv),length(woe_conv));
        % Reward on Tin
        pick = find(TinRew & correct & num_accum==gi);
        pick_cor = pick;
        if isempty(pick)
            TinRewCor{gi} = zero_temp;
        else
            cor_woe = cumLLR{1}(pick,gi);
            cor_woe_all = [cor_woe_all;cor_woe];
            cor_fr = frMtrx{1}(pick,gi);
            cor_fr_all = [cor_fr_all;cor_fr];
            TinRewCor{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        
        pick = find(TinRew & ~correct & num_accum==gi);
        pick_res = pick;
        if isempty(pick)
            TinRewErr{gi} = zero_temp;
        else
            TinRewErr{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        
        pick = find(TinRew & num_accum~=gi);
        if isempty(pick)
            TinRewRes{gi} = zero_temp;
        else
            res_woe = cumLLR{1}(pick,gi);
            res_woe_all = [res_woe_all;res_woe];
            res_fr = frMtrx{1}(pick,gi);
            res_fr_all = [res_fr_all;res_fr];
            
            res_woe = cumLLR{1}(pick,gi);
            res_fr = frMtrx{1}(pick,gi);
            TinRewRes{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        if ~isempty(pick_cor) && ~isempty(pick_res)
            [AUC hit fa] = calcROC(cor_woe, res_woe);
            AUC_woe(gi) = AUC;
            [AUC hit fa] = calcROC(cor_fr, res_fr);
            AUC_fr(gi) = AUC;
        end
        
        
        % Reward on Tout
        pick = find(ToutRew & correct & num_accum==gi);
        if isempty(pick)
            ToutRewCor{gi} = zero_temp;
        else
            ToutRewCor{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        
        pick = find(ToutRew & ~correct & num_accum==gi);
        if isempty(pick)
            ToutRewErr{gi} = zero_temp;
        else
            ToutRewErr{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        
        pick = find(ToutRew & num_accum~=gi);
        if isempty(pick)
            ToutRewRes{gi} = zero_temp;
        else
            ToutRewRes{gi} = hist2(cumLLR{1}(pick,gi),frMtrx{1}(pick,gi),woe_conv,fr_conv);
        end
        figure(fig)
        subplot(3,10,gi)
        imagesc(woe_conv/10,fr_conv,TinRewCor{gi})
        axis xy
        xlim([-4 4]);
        subplot(3,10,10+gi)
        imagesc(woe_conv/10,fr_conv,TinRewRes{gi})
        axis xy
        xlim([-4 4]);
        subplot(3,10,20+gi)
        imagesc(woe_conv/10,fr_conv,TinRewErr{gi})
        axis xy
        xlim([-4 4]);
        
        figure(fig+1)
        subplot(3,10,gi)
        imagesc(woe_conv/10,fr_conv,fliplr(ToutRewErr{gi}))
        axis xy
        xlim([-4 4]);
        subplot(3,10,10+gi)
        imagesc(woe_conv/10,fr_conv,fliplr(ToutRewRes{gi}))
        axis xy
        xlim([-4 4]);
        subplot(3,10,20+gi)
        imagesc(woe_conv/10,fr_conv,fliplr(ToutRewCor{gi}))
        axis xy
        xlim([-4 4]);
        
        figure(fig+2)
        subplot(2,10,gi);hold on
        plot(woe_conv/10,sum(TinRewCor{gi},1)./sum(sum(TinRewCor{gi},1)),'-r')
        plot(woe_conv/10,sum(TinRewRes{gi},1)./sum(sum(TinRewRes{gi},1)),'-k')
        xlim([-4 4]);
        subplot(2,10,10+gi);hold on
        % fr_axis = 1:1:101;
        plot(fr_conv,sum(TinRewCor{gi},2)./sum(sum(TinRewCor{gi},2)),'-r')
        plot(fr_conv,sum(TinRewRes{gi},2)./sum(sum(TinRewRes{gi},2)),'-k')
        xlim([0 100])
    end
    
    [AUC_woe_all hit fa] = calcROC(cor_woe_all(~isnan(cor_woe_all)), res_woe_all(~isnan(res_woe_all)));
    [AUC_fr_all hit fa] = calcROC(cor_fr_all(~isnan(cor_fr_all)), res_fr_all(~isnan(res_fr_all)));
    info.AUC_woe = AUC_woe_all;
    info.AUC_fr = AUC_fr_all;
    
    figure(fig+3); hold on
    plot(1:10,AUC_woe,'-bo');
    plot(1:10,AUC_fr,'-ro');
    
    FH4 = figure(fig+4);
    set(FH4,'position',[100 100 300 100]);
    
    
    num_fig = 4;
    fig = fig+num_fig;
end

%% Variance analysis
if fig_switch(21)
    num_fig = 1;
    figure(fig);clf;hold on
    uniqueNum = 4:8; %unique(num_accum);
    fr_var = zeros(length(uniqueNum),uniqueNum(end));
    % creating color index
    map = colormap;
    interval = floor((64-1)/(length(uniqueNum)-1));
    plot_color = map(1:interval:(1+interval*(length(uniqueNum)-1)),:);
    for i = 1:length(uniqueNum)
        for gi = 1:2
            pick = logical(num_accum==uniqueNum(i) & choice == gi);
            for k = 1:uniqueNum(i)+2
                fr_var(i,gi,k)=var(frMtrx{1}(pick,k));
            end
        end
        subplot(2,6,i)
        plot(1:i+2,squeeze(fr_var(i,1,1:i+2)),'color',plot_color(i,:))
        plot(1:i+2,squeeze(fr_var(i,1,1:i+2)),'color',plot_color(i,:))
    end
    xlabel('Number of shapes','FontSize',15,'FontWeight','bold'); 
    ylabel('Variance of FR','FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',15,'FontWeight','bold','Box','OFF','TickDir','out')
    hold off
    fig = fig+num_fig;
end
%%

return