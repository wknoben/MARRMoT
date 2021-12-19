load all_models.mat
all_outputs = cell(numel(all_models),1);
all_times = zeros(numel(all_models),1);
for i = 1:numel(all_models)
    if isempty(all_outputs{i})
        tic
        all_outputs{i} = run_with_rand_par(all_models(i));
        all_times(i) = toc;
    end
end