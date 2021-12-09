function [example_ind] = get_indices_at_y_level(targets, y)

for i = 1:length(targets)
    t = targets(i);
    abs_diff = abs(y - t);
    [~, i_min] = min(abs_diff);
    example_ind(i) = i_min;
end

end

