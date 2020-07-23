% Import model results
radi = load('../results/IC_W29.mat');

% Get names of modelled variables
vars = fieldnames(radi);
vars = vars(~ismember(vars, {'savetimes' 'depths'}));

% Define timestep colours
clr = turbo(numel(radi.savetimes));

% Plot all results
figure(1)
clf
for v = 1:numel(vars)
  subplot(4, 5, v)
  hold on
  for t = 1:numel(radi.savetimes)
    plot(radi.depths, radi.(vars{v})(:, t), 'color', clr(t, :))
  end  % for t
  xlabel('Depth / cm')
  ylabel(vars{v})
  grid on
end  % for v
