cd('/Users/cameronrichardson/Documents/School/Thesis/Code');

% Read the file using readtable with a comma delimiter
data = readtable('ID_pK1.txt', 'Delimiter', ',');

% Assuming 'ID' is a column from your 'data' table
ID = data(:, 1);  % Column containing ID strings
pK1 = data{:, 2}; % Convert to numeric array (if it's not already)

% Convert the ID table column to a cell array for string manipulation
ID = table2cell(ID); 

% Initialize the group cell array
group = cell(size(ID, 1), 1);  % Use size() for arrays

% Grouping based on patterns in ID strings
for i = 1:size(ID, 1)
    if contains(ID{i}, 'BP2')
        group{i} = 'BP2';
    elseif contains(ID{i}, 'BP3')
        group{i} = 'BP3';
    elseif contains(ID{i}, 'BP4')
        group{i} = 'BP4';
    elseif contains(ID{i}, 'BP5')
        group{i} = 'BP5';
    else
        group{i} = 'pK1';
    end
end

% delta pK1 loop
delta = zeros(size(pK1));  % Pre-allocate numeric array for delta values

for i = 1:size(pK1, 1)
    delta(i) = pK1(1) - pK1(i);  % Difference with first element in pK1
end

% Convert group to categorical for easy handling in plots
group = categorical(group);

% Calculate the mean and standard deviation of pK1
pk1_true = pK1(2:16);
mean_pK1 = mean(pk1_true);
std_pK1 = std(pk1_true);





% Calculate the confidence intervals (95%)
n = length(pk1_true);  % Number of samples
conf_level = 0.95;  % Confidence level
alpha = 1 - conf_level;  % Significance level
t_score = tinv(1 - alpha/2, n - 1);  % t-score for 95% confidence interval

margin_of_error = t_score * (std_pK1 / sqrt(n));  % Margin of error
CI_lower = mean_pK1 - margin_of_error;  % Lower bound of confidence interval
CI_upper = mean_pK1 + margin_of_error;  % Upper bound of confidence interval

% Create a simple plot with filled markers based on groups
figure;
hold on;  % Hold on to plot multiple groups

% Create a color map
colors = lines(numel(categories(group)));  % Get distinct colors for each group

% Plot filled markers for each group
for g = categories(group)'
    idx = group == g;  % Get index for the current group
    plot(pK1(idx), delta(idx), 'o', 'Color', colors(find(strcmp(categories(group), g)), :), ...
        'MarkerFaceColor', colors(find(strcmp(categories(group), g)), :), ... % Fill color
        'MarkerSize', 8, 'DisplayName', char(g));
end

% Plot mean and confidence intervals
plot([mean_pK1, mean_pK1], [min(delta), max(delta)], 'k--', 'LineWidth', 2, 'DisplayName', 'Mean pK1');  % Mean line
plot([CI_lower, CI_lower], [min(delta), max(delta)], 'r--', 'LineWidth', 1.5, 'DisplayName', '95% CI Lower');  % Lower CI
plot([CI_upper, CI_upper], [min(delta), max(delta)], 'r--', 'LineWidth', 1.5, 'DisplayName', '95% CI Upper');  % Upper CI

% Add labels and title
xlabel('pK1');
ylabel('Delta pK1');
title('Filled Markers Plot of pK1 vs Delta pK1 with Confidence Intervals');
legend('show');  % Show legend for the groups
hold off;  % Release the hold





% Remove data points that fall outside the confidence intervals
valid_indices = pk1_true >= CI_lower & pk1_true <= CI_upper;
filtered_pK1 = pk1_true(valid_indices);
filtered_delta = delta(valid_indices);
filtered_group = group(valid_indices);

% Recalculate mean and standard deviation for the filtered data
new_mean_pK1 = mean(filtered_pK1);
new_std_pK1 = std(filtered_pK1);


pK1_BP2mean = mean(pK1(2:4));
pK1_BP2std = std(pK1(2:4));

pK1_BP3mean = mean(pK1(6:10));
pK1_BP3std = std(pK1(6:10));



pK1_BP5mean = mean(pK1(11:16));
pK1_BP5std = std(pK1(11:16));
% Recalculate confidence intervals for filtered data
n_filtered = length(filtered_pK1);  % Number of samples in filtered data
if n_filtered > 1  % Ensure there are enough points for t-test
    new_t_score = tinv(1 - alpha/2, n_filtered - 1);  % t-score for 95% confidence interval
    new_margin_of_error = new_t_score * (new_std_pK1 / sqrt(n_filtered));  % New margin of error
    new_CI_lower = new_mean_pK1 - new_margin_of_error;  % New lower CI
    new_CI_upper = new_mean_pK1 + new_margin_of_error;  % New upper CI
else
    new_CI_lower = new_mean_pK1;  % If only one point, CI is the mean
    new_CI_upper = new_mean_pK1;
end

% Create a simple plot with filled markers based on groups
figure;
hold on;  % Hold on to plot multiple groups

% Create a color map
colors = lines(numel(categories(filtered_group)));  % Get distinct colors for each group

% Plot filled markers for each group
for g = categories(filtered_group)'
    idx = filtered_group == g;  % Get index for the current group
    plot(filtered_pK1(idx), filtered_delta(idx), 'o', 'Color', colors(find(strcmp(categories(filtered_group), g)), :), ...
        'MarkerFaceColor', colors(find(strcmp(categories(filtered_group), g)), :), ... % Fill color
        'MarkerSize', 8, 'DisplayName', char(g));
end

% Plot new mean and confidence intervals
plot([new_mean_pK1, new_mean_pK1], [min(filtered_delta), max(filtered_delta)], 'k--', 'LineWidth', 2, 'DisplayName', 'Mean pK1');  % Mean line
plot([new_CI_lower, new_CI_lower], [min(filtered_delta), max(filtered_delta)], 'r--', 'LineWidth', 1.5, 'DisplayName', '95% CI Lower');  % Lower CI
plot([new_CI_upper, new_CI_upper], [min(filtered_delta), max(filtered_delta)], 'r--', 'LineWidth', 1.5, 'DisplayName', '95% CI Upper');  % Upper CI

% Add labels and title
xlabel('pK1');
ylabel('Delta pK1');
title('Filtered Filled Markers Plot of pK1 vs Delta pK1 with Confidence Intervals');
legend('show');  % Show legend for the groups
hold off;  % Release the hold