clear; clc;

% read in protist count data
[num words] = xlsread('Halteria_PB_5-6-22.xlsx','counts');

% read in plaque assay data
[num2 words2] = xlsread('Halteria_PB_5-6-22.xlsx','Plaque assay');

%% this section plots cell counts and PFUs

days = [0,1,2];
number_0 = num(:,6);
number_1 = num(:,10);
number_2 = num(:,13);
counts = [number_0, number_1, number_2];
treatment = words(2:25,3);
species = words(2:25,2);

% pull out and plot halteria data
HA_rows_V = find(ismember(treatment,'V') & ismember(species,'Halteria'));
HA_rows_NV = find(ismember(treatment,'NV') & ismember(species,'Halteria'));

figure(1); clf(1);
subplot(2,2,3);
box on; hold on;
    for i = 1:length(HA_rows_V)
        plot(days,counts(HA_rows_V(i),:),'-','Color',[0.5 0.5 0.5]);
    end
        h1 = plot(days,mean(counts(HA_rows_V,:)),'-k','LineWidth',2);
    for i = 1:length(HA_rows_NV)
        plot(days,counts(HA_rows_NV(i),:),'-','Color',[0.5 0.5 1]);
    end
        h2 = plot(days,mean(counts(HA_rows_NV,:)),'-b','LineWidth',2);
    ylim([0 30]);
    xticks([0,1,2]);
    ylabel('Number of cells','FontSize',12);
    xlabel('Days','FontSize',12);
    text(0.1,27,'C','FontSize',12);
    legend([h1 h2],'With virus','Without virus','Location','Best');

% pull out and plot paramecium data
PB_rows_V = find(ismember(treatment,'V') & ismember(species,'Paramecium bursaria'));
PB_rows_NV = find(ismember(treatment,'NV') & ismember(species,'Paramecium bursaria'));

subplot(2,2,4);
box on; hold on;
    for i = 1:length(PB_rows_V)
        plot(days,counts(PB_rows_V(i),:),'-','Color',[0.5 0.5 0.5]);
    end
        plot(days,mean(counts(PB_rows_V,:)),'-k','LineWidth',2);
    for i = 1:length(PB_rows_NV)
        plot(days,counts(PB_rows_NV(i),:),'-','Color',[0.5 0.5 1]);
    end
        plot(days,mean(counts(PB_rows_NV,:)),'-b','LineWidth',2);
    xticks([0,1,2]);
    text(0.1,27,'D','FontSize',12);
    ylim([0 30]);
    xlabel('Days','FontSize',12);

% now do virus counts
PFUs = num2(:,3);

subplot(2,2,1);
box on;
    for i = 1:6
        semilogy([0 2],[1.11e7 PFUs(i)],'-','Color',[0.5 0.5 0.5]); hold on;
    end
    ylim([1e5 2e7]);
    ylabel('PFUs mL^{-1}','FontSize',12);
    text(0.1,2e5,'A','FontSize',12);
    xticks([0,1,2]);
    title('\itHalteria');

subplot(2,2,2);
box on;
    for i = 7:12
        semilogy([0 2],[1.11e7 PFUs(i)],'-','Color',[0.5 0.5 0.5]); hold on;
    end
    ylim([1e5 2e7]);
    title('\itParamecium');
    text(0.1,2e5,'B','FontSize',12);
    xticks([0,1,2]);

    shg;


%% read in fitted ODE parameters
[num3 words3] = xlsread('Bootstrapped_parameters.xls','Sheet1');

scr = num3(:,2);
ce = num3(:,3);
hstart = num3(:,4);

HA_counts = counts(HA_rows_V,:);
HA_PFUs = PFUs(1:6)./0.3;

% ODE setup
tspan = 0:0.1:2; % start and end times

figure(3);clf(3);
subplot(1,2,1);
hold on; box on;
    jitters = 1 + 0.02*randn(size(HA_counts,1),size(HA_counts,2));
    plot(jitters.*days,jitters.*HA_counts,'ok','MarkerFaceColor',[0.5 0.5 0.5]);
    xlim([-0.1 2.1]);
    text(0.1,27,'A','FontSize',12);
    xlabel('Days','FontSize',12);
    ylabel('Cells mL^{-1}','FontSize',12);
    title('Halteria');

subplot(1,2,2);
box on;
    semilogy([2 2 2 2 2 2],HA_PFUs,'ok','MarkerFaceColor',[0.5 0.5 0.5]); hold on; 
    semilogy([0 0 0 0 0 0],3.7e7,'ok','MarkerFaceColor',[0.5 0.5 0.5]); 
    ylim([1e5 4e7]);    
    xlabel('Days','FontSize',12);
    ylabel('Virus particles mL^{-1}','FontSize',12);
    title('Chlorovirus');
    text(0.1,2e5,'B','FontSize',12);
    xticks([0,1,2]);
    xlim([-0.1 2.1]);    
    
    for i = 1:length(scr)        
        y0 = [3.7e7, hstart(i)];
        ode = @(t,y) Virovore_virus_model(t,y,scr(i),ce(i));
        [t1,y1] = ode45(ode, tspan, y0); % return time and population density vectors
        figure(3);
        subplot(1,2,1);
            plot(t1,y1(:,2),'-');
        subplot(1,2,2);
            plot(t1,y1(:,1),'-');            
    end

    shg;
