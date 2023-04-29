function [time, injs, y_conc, x_conc] = get_data(true_params, config)
    %% Create a two compartment model structure
    pk_construct = PKModelDesign;
    
    addCompartment(pk_construct, 'Central',...
                   'DosingType', 'Infusion', ...
                   'EliminationType', 'linear', ...
                   'HasResponseVariable', true);
    
    addCompartment(pk_construct, 'Peripheral',...
                   'DosingType', '', ...
                   'HasResponseVariable', true);
    
    [modelObj, pkmap] = pk_construct.construct;
    
    %% Set Parameter Values
    % Set compartment volumes
    modelObj.Compartments(1).Capacity = true_params.Central;
    modelObj.Compartments(2).Capacity = true_params.Peripheral;
    
    % Set important parameters reasonable values and others 0 and convert units
    important_params = pkmap.Estimated;
    for i = 1:length(modelObj.Parameters)
        current_unit = get(modelObj.Parameters(i), 'Units');
        new_unit = strrep(current_unit, 'hour', 'minute');
        set(modelObj.Parameters(i), 'Units', new_unit);
        is_important = cellfun(@(x) strcmp(x,modelObj.Parameters(i).Name), important_params);
        if ~any(is_important)
            set(modelObj.Parameters(i), 'Value', 0);
        else
            name_matched = important_params{is_important};
            set(modelObj.Parameters(i), 'Value', true_params.(name_matched));
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
    configset.StopTime = 400;
    
    %% Set Dose Specs
    doseObj = adddose(modelObj,'InfusionDose');
    doseObj.TargetName = 'Drug_Central';
    doseObj.StartTime = 20;
    doseObj.Amount = 10;
    doseObj.Rate = 0.5;
    doseObj.Interval = 120;
    doseObj.RepeatCount = config.repeat_count-1;
    doseObj.AmountUnits = 'micromole';
    doseObj.TimeUnits = 'minute';
    doseObj.RateUnits = 'micromole/minute';
    
    %% Simulate model and add noise
    % Simulate and create injections vector
    simdata = sbiosimulate(modelObj, doseObj);
    if config.uniform_sampling
        Tsample = config.sampling_time;
        time = (0:Tsample:configset.StopTime)';
        newSimData = resample(simdata,time, 'spline');
        x_conc = newSimData.Data;
    else
        time = simdata.Time;
        x_conc = simdata.Data;
    end
    
    injection_length = doseObj.Amount/doseObj.Rate;
    injs = zeros(size(time));
    start_times = doseObj.StartTime+(0:doseObj.RepeatCount)*doseObj.Interval;
    for i=1:length(start_times)
        nnz_inj = time>=start_times(i) &...
                  time<=start_times(i)+injection_length;
        injs(nnz_inj) = doseObj.Rate;
    end
    
    % Add noise of particular SNR
    SNR = 30;
    rng(1)
    noise = normrnd(0, 1, size(x_conc));
    
    sig_pow = diag(x_conc'*x_conc)/size(x_conc, 1);
    y_conc = x_conc+noise*sqrt(diag(sig_pow/SNR));
end