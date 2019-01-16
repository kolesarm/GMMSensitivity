function figures
    addpath(genpath('../external'))
    load ../output/transp_identification_blp
    load ../external/data/bootstrap_std_sen

    fid = fopen('param_crosswalk.csv');
    fgetl(fid);
    param_crosswalk = textscan(fid,'%s %s','delimiter',',');
    fclose(fid);

    fid = fopen('moment_crosswalk.csv');
    fgetl(fid);
    mom_crosswalk = textscan(fid,'%s %s','delimiter',',');
    fclose(fid);

    opts = VisualizeTableOptions('pad_width', 4, 'pad_height', 2);
    opts.label_props.row.HorizontalAlignment = 'right';
    opts = opts.shift_label('row', 0.45, 0);
    opts = opts.shift_label('col', -.25, .4);

    % Display sensitivity of blp1995 parameters to blp1995 moments
    rows = get_crosswalk_names(blp1995.paramlist, param_crosswalk);
    rows = strrep(rows, 'Air conditioning', 'AC');
    rows = strrep(rows, 'random coefficient', 'standard deviations');

    demand_moments = strncmp(blp1995.momlist, 'demand', 6);
    cols_demand = get_crosswalk_names(blp1995.momlist(demand_moments), mom_crosswalk);
    cols_demand = strcat(cols_demand, ' \times \xi');
    table = blp1995.standardized_sensitivity(:, demand_moments);
    create_table(table, cols_demand, rows, 'blp1995_demand', opts);

    se_opts = opts;
    se_opts.se_mat = bootstrap_se_std_sen(:, demand_moments);
    create_table(table, cols_demand, rows, 'blp1995_demand_with_se', se_opts);

    supply_moments = strncmp(blp1995.momlist, 'supply', 6);
    cols_supply = get_crosswalk_names(blp1995.momlist(supply_moments), mom_crosswalk);
    cols_supply = strcat(cols_supply, ' \times \omega');
    table = blp1995.standardized_sensitivity(:, supply_moments);
    create_table(table, cols_supply, rows, 'blp1995_supply', opts);

    se_opts.se_mat = bootstrap_se_std_sen(:, supply_moments);
    create_table(table, cols_supply, rows, 'blp1995_supply_with_se', se_opts);

    % Sample sensitivity tables
    table = blp1995.standardized_sample_sensitivity(:, demand_moments);
    create_table(table, cols_demand, rows, 'blp1995_demand_sample', opts);

    table = blp1995.standardized_sample_sensitivity(:, supply_moments);
    create_table(table, cols_supply, rows, 'blp1995_supply_sample', opts);

    exit
end

function rows = get_crosswalk_names(varnames, crosswalk, sufficiency)
    rows = cell(length(varnames),1);
    for i = 1:length(varnames)
        j = strncmp(varnames(i), crosswalk{1}, 20);
        name = crosswalk{2}(j);
        if nargin > 2
            rows{i} = sprintf('%s\n(sufficiency = %.3f)', name{:}, sufficiency(i));
        else
            rows{i} = name{:};
        end
    end
end

function create_table(table, cols, rows, name, opts)
    opts.precision = 2;
    opts.cell_width = 0.5;
    opts.cell_height = 0.4;
    opts = opts.shift_label('row', 0.35, 0);
    opts = opts.shift_label('col', -.2, .2);
    opts.label_props.col.Rotation = 25;

    % Display table of sensitivity values
    visualize_table(table, rows, cols, opts);

    filename = sprintf('../temp/%s_table.png', name);
    saveas(gcf, filename, 'png');

    filename = sprintf('../output/%s_table.eps', name);
    print(filename, '-depsc2');
end

function create_heatmap(table, cols, rows, name, opts)
    % Make color map using reversed data (so that large values are darker)
    hm_data = -table;
    opts.binary_color = false;
    opts.colormap_scale = 'bone';
    opts.gridlines = false;
    opts.tabletext_props.Visible = 'off';

    visualize_table(hm_data, rows, cols, opts);

    title('Standardized sensitivity heat map', 'FontSize', 14, ...
                                               'FontWeight', 'bold', ...
                                               'Position', [3.5, 0]);

    % Legend for heat map
    cbar = colorbar();
    neg_ticks_char = get(cbar, 'YTickLabel');
    ticks = zeros(1, length(neg_ticks_char));
    for i = 1:length(neg_ticks_char)
        ticks(i) = -str2double(neg_ticks_char(i, :));
    end
    colorbar('location', 'eastoutside', 'YTickLabel', ticks);

    filename = sprintf('../temp/%s_heatmap.png', name);
    saveas(gcf, filename, 'png');

    filename = sprintf('../output/%s_heatmap.eps', name);
    print(filename, '-depsc2');
end
