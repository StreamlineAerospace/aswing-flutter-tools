clc; 
clear variables; 
close all;
warning off;

% Script creates V-F and V-G plots from ASWING root-locus output
% Choose input file to read eigenvalues
% Choose file with airspeed point list

rtfilename = 'Ptest.e00'; % root locus data filename
speedfilename = 'Ptest.t'; % airspeed filename

%*****************************************************************

data = readtable(rtfilename,"FileType","text");
airspeeds = readtable(speedfilename,"FileType","text",'NumHeaderLines',2);
airspeeds.VIAS_ref = airspeeds.V;
airspeeds.V = airspeeds.t;

eigenvalues = floor(height(data)/length(airspeeds.VIAS_ref))+3;
point_index(1) = 1;
for i = 1:height(data)
    if length(point_index) < data.x_(i)
        point_index(data.x_(i)) = i;
    end
end
point_index(end+1) = height(data)+1;

% Read Data
first_flutter=nan(4,1);
j=1;
for i = 1:height(data)
    f(j,i+1-point_index(j)) = table2array(data(i,3));
    sigma(j,i+1-point_index(j)) = table2array(data(i,2));
    theta_temp = atan(table2array(data(i,3))/table2array(data(i,2)));
    if theta_temp < 0
        theta_temp = pi() + theta_temp;
    end
    G(j,i+1-point_index(j)) = 2*cos(theta_temp);
    if isnan(first_flutter(1)) && sigma(j,i+1-point_index(j)) > 0
        v_index = double(data.x_(find(ismember(data.Eigenvalue,f(j,i+1-point_index(j))))));
        first_flutter = [sigma(j,i+1-point_index(j)),f(j,i+1-point_index(j)),airspeeds.VIAS_ref(v_index),airspeeds.V(v_index)];
    end
    if point_index(j+1) == i+1
        j = j+1;
    end     
end

% Fix Frequencies and G values
f(f == 0) = NaN;
G(G == 0) = NaN;
for i = 1:length(f)
    tempf = f(i,:);
    tempsigma = sigma(i,:);
    tempg = G(i,:);
    [~,sort_index] = sort(f(i,:),2);
    for j = 1:length(tempf)
        f(i,j) = tempf(sort_index(j));
        sigma(i,j) = tempsigma(sort_index(j));
        G(i,j) = tempg(sort_index(j));
    end
end


% for i = 2:length(f)
%     new_index = zeros(length(f(i,:)),1);
%     k=0;
%     for j = 1:length(f(i,:))
%         if isnan(f(i,j)) == 0
%             [minval,temp] = min(abs(f(i,j)-f(i-1,:)));  
% %             if ismember(temp,new_index) == 1
% %                 k = k+1;
% %             end
%             new_index(j+k) = temp;
%         end
%     end
%     newf = nan(length(f(i,:)),1);
%     newsigma = nan(length(f(i,:)),1);
%     for j = 1:length(f(i,:))
%         if new_index(j) ~= 0
%             newf(new_index(j)) = f(i,j);
%             newsigma(new_index(j)) = sigma(i,j);
%         end
%     end
%     f(i,:) = newf;
%     sigma(i,:) = newsigma;
% 
%     new_index = zeros(length(G(i,:)),1);
%     k=0;
%     for j = 1:length(G(i,:))
%         if isnan(G(i,j)) == 0
%             [minval,temp] = min(abs(G(i,j)-G(i-1,:)));  
%             if ismember(temp,new_index) == 1
%                 k = k-1;
%             end
%             new_index(j+k) = temp;
%         end
%     end
%     newG = nan(length(G(i,:)),1);
%     for j = 1:length(G(i,:))
%         if new_index(j) ~= 0
%             newG(new_index(j)) = G(i,j);
%         end
%     end
%     G(i,:) = newG;
% end


% Find first flutter crossing
first_flutter=nan(4,1);
for i = 1:width(f)
    for j = 2:length(f(:,i))
        if isnan(first_flutter(1)) && 0>sigma(j-1,i)*sigma(j,i) && abs(f(j-1,i)-f(j,i))<0.5
            v_index = double(data.x_(find(ismember(data.Eigenvalue,f(j,i)))));
            first_flutter = [sigma(j,i),f(j,i)/(2*pi()),airspeeds.V(v_index),airspeeds.VIAS_ref(v_index)];
        end
    end
end

fprintf('First flutter occurs at sigma = ' + string(first_flutter(1)) + ' [Hz] and ω = ' + string(first_flutter(2)));
fprintf('\nThis occurs at an indicated airspeed of ' + string(first_flutter(3)) + ' m/s.\n');

% Plot Flutter Results
figure();
sgtitle('Root Locus Plot');
hold on;
for i = 1:length(f)
    plot(sigma(i,1:point_index(i+1)-point_index(i)),f(i,1:point_index(i+1)-point_index(i))/(2*pi()),'.','MarkerSize',10);
end
if first_flutter(2) >0
    plot(first_flutter(1),first_flutter(2),'*r','MarkerSize',15);
end
hold off;
line([0,0], ylim,'LineStyle','--', 'Color', 'k', 'LineWidth', 1); % Draw line for X axis.
grid on;
axis on;
xlabel('σ [Hz]');
ylabel('Frequency [Hz]');

% V-F plot
figure();
sgtitle('Propsperity 1 V-F and V-G Plots');
subplot(2,1,1);           
% title('V-F Plot');
hold on;
for i = 1:width(f)
    end_range(i) = length(nonzeros(f(:,i)));
    plot(airspeeds.VIAS_ref(1:end_range(i)),f(1:end_range(i),i)/(2*pi()),'.','MarkerSize', 10);
    legend_string(i) = 'Eigenvalue #' + string(i);
end
hold off;
grid on;
set(gca,'XAxisLocation','top');
xlim([min(airspeeds.VIAS_ref)-10 max(airspeeds.VIAS_ref)+10]);
xlabel('V-F Plot: Indicated Airspeed (m/s)');
ylabel('Frequency [Hz]');
legend(legend_string,'Location','southeast');

% V-G plot
subplot(2,1,2);
% title('V-G Plot');
xlabel('True Airspeed (m/s)');
xlim([min(airspeeds.V)-10 max(airspeeds.V)]+10);
legend(legend_string,'Location','southeast');
hAx(1)=gca;
set(gca,'ytick',[],'XColor','b');
grid on;
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','left');
hold(hAx(2),'on')
xlim([min(airspeeds.VIAS_ref)-10 max(airspeeds.VIAS_ref)+10]);
ylim([-0.7 0.05]);
xlabel('V-G Plot: Indicated Airspeed (m/s)');
ylabel('Damping');
hold on;
for i = 1:width(G)
    plot(airspeeds.VIAS_ref(1:end_range(i)),G(1:end_range(i),i),'.','MarkerSize', 10);
end
hold off;
grid on;
line(xlim, [0,0],'LineStyle','--', 'Color', 'k', 'LineWidth', 1); % Draw line for X axis.
legend(legend_string,'Location','southwest');


