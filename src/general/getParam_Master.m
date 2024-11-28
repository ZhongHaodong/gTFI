function value = getParam_Master(params, name, defaultValue, master_mode)
    if master_mode
        value = getOrDefault(params, name, defaultValue);
    else
        value = defaultValue;
    end
end
