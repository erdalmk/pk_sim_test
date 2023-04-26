%% Create a two compartment model structure
pk_construct = PKModelDesign;

comp1 = addCompartment(pk_construct, 'Central',...
                       'DosingType', 'Infusion', ...
                       'EliminationType', 'linear', ...
                       'HasResponseVariable', true);

comp2 = addCompartment(pk_construct, 'Peripheral',...
                       'DosingType', '', ...
                       'HasResponseVariable', true);

[modelObj, pkmap] = pk_construct.construct;

%% Set Parameter Values
% Set compartment volumes
modelObj.Compartments(1).Capacity = 5e-2;
modelObj.Compartments(2).Capacity = 3e-2;

% Set important parameters reasonable values and others 0 and convert units
important_params = pkmap.Estimated;
important_param_vals = [modelObj.Compartments(1).Capacity,...
                        3e-2,...
                        modelObj.Compartments(2).Capacity,...
                        3e-4];
for i = 1:length(modelObj.Parameters)
    current_unit = get(modelObj.Parameters(i), 'Units');
    new_unit = strrep(current_unit, 'hour', 'minute');
    set(modelObj.Parameters(i), 'Units', new_unit);
    is_important = cellfun(@(x) strcmp(x,modelObj.Parameters(i).Name), important_params);
    if ~any(is_important)
        set(modelObj.Parameters(i), 'Value', 0);
    else
        set(modelObj.Parameters(i), 'Value', important_param_vals(is_important));
    end
end

for i = 1:length(modelObj.Species)
    current_unit = get(modelObj.Species(i), 'Units');
    new_unit = strrep(current_unit, 'milligram', 'micromole');
    set(modelObj.Species(i), 'Units', new_unit);
end

%% Set model configurations for simulation
configset = getconfigset(modelObj);
configset.CompileOptions.UnitConversion = true;
configset.TimeUnits = 'minute';
configset.StopTime = 600;

%% Set Dose Specs
doseObj = adddose(modelObj,'InfusionDose');
doseObj.TargetName = 'Drug_Central';
doseObj.StartTime = 20;
doseObj.Amount = 10;
doseObj.Rate = 1;
doseObj.AmountUnits = 'micromole';
doseObj.TimeUnits = 'minute';
doseObj.RateUnits = 'micromole/minute';

%% Simulate model and add noise
% Simulate
[time, x_conc,names] = sbiosimulate(modelObj, doseObj);

% Add noise of particular SNR
SNR = 50;
rng(1)
noise = normrnd(0, 1, size(x_conc));

sig_pow = sqrt(diag(x_conc'*x_conc));
y_conc = x_conc+noise*diag(sig_pow/SNR);


%% Plot simulated response
figure
plot(time, y_conc, '.')
hold on
set(gca,'ColorOrderIndex',1)
plot(time, x_conc)
xlim(time([1, end]))
ylim([-20, 200])
set(gca,'TickDir', 'out', 'box', 'off')