function value = getOrDefault(struct, field, defaultValue)
    if isfield(struct, field)
        value = struct.(field);
    else
        value = defaultValue;
    end
end
