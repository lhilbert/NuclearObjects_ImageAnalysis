classdef DatasetFileManager
    %Handles data files in hierarchical folder structures
    %
    %   The DATASETFILEMANAGER class provides methods to work with data files
    %   from experiment campaigns, that are collected in hierarchical folder
    %   structures, with experiment parameters encoded arbitrarily in the folder
    %   and file names.
    %
    %   Start with defining an arbitrary number of experiment parameters (using
    %   addParameter()) and the specific expressions available for every
    %   parameter (using addLabel()). For each parameter expression, the class
    %   stores a unique label, a filepath component that specifies the
    %   expression in the original data file paths, and optionally a filepath
    %   component that specifies the expression in a set of extracted data
    %   files.
    %
    %   The OriginalFilepathPattern and ExtractedFilepathPattern properties
    %   specify the filepath patterns, in which the parameter-specifying
    %   filepath components are provided using the placeholders '{P1}', '{P2}',
    %   '{P3}', ..., the number corresponding with the order in which parameters
    %   were added.

    properties
        OriginalFilepathPattern  (1,1) string = "";
        ExtractedFilepathPattern (1,1) string = "";
    end

    properties (Access=private)
        ParameterTables (:,1) cell = cell(0,1);
        ParameterNames  (:,1) string = strings(0,1);
    end

    methods
        function obj = DatasetFileManager()
            %Construct new DATASETFILEMANAGER instance
        end

        function tab = getParameter(obj, par_name)
            %Retrieve parameter table
            par_ind = obj.findParameter(par_name);
            if isempty(par_ind)
                tab = [];
            else
                tab = obj.ParameterTables{par_ind};
            end
        end

        function obj = addParameter(obj, par_name)
            %Add new parameter table
            assert(isempty(obj.findParameter(par_name)), ...
                "Parameter '%s' already exists. Use unique parameter names.", par_name);
            obj.ParameterNames(end+1) = par_name;
            obj.ParameterTables(end+1) = {DatasetFileManager.initializeParameterTable(par_name)};
        end

        function tab = getLabel(obj, par_name, label)
            %Retrieve parameter expression row from parameter table
            [par_ind, label_ind] = obj.findlabel(par_name, label);
            if isempty(par_ind) || isempty(label_ind)
                tab = [];
            else
                tab = obj.ParameterTables{par_ind}(label_ind,:);
            end
        end

        function obj = addLabel(obj, par_name, label, org_filepath_comp, ext_filepath_comp)
            %Add new parameter expression to parameter table
            par_ind = obj.findParameter(par_name);
            assert(~isempty(par_ind), "Parameter '%s' not found. Add it first with addParameter().", par_name);
            assert(isempty(find(obj.ParameterTables{par_ind}{:,1} == label, 1)), ...
                compose("Parameter '%s' already contains label '%s'. Use unique labels.", par_name, label));
            obj.ParameterTables{par_ind}(end+1,:) = {label, org_filepath_comp, ext_filepath_comp};
        end

        function tab = compileParameterCombinations(obj)
            %Create all parameter expression combinations
            combos = obj.labelIndexCombos();
            tab = obj.initializeFileTable();
            tab = removevars(tab, 1);
            for p = 1:height(combos)
                tab_values = obj.retrieveParameterValues(combos{p,:}, 1:3);
                tab{end+1, :} = tab_values; %#ok<AGROW>
            end
        end

        function tab = compileOriginalFileTable(obj)
            %Retrieve full file list matching OriginalFilepathPattern
            assert(~isempty(obj.OriginalFilepathPattern), ...
                "OriginalFilepathPattern not set. Assign it before compiling file table.");
            n_pars = length(obj.ParameterTables);
            combos = obj.labelIndexCombos();
            tab = obj.initializeFileTable();
            filepath_placeholders = compose("{P%d}", 1:n_pars);
            for p = 1:height(combos)
                filepath_values = obj.retrieveParameterValues(combos{p,:}, 2);
                file_pattern = replace(obj.OriginalFilepathPattern, filepath_placeholders, filepath_values);
                files = dir(file_pattern);
                files = fullfile(string({files.folder}), string({files.name}))';
                n_files = length(files);
                tab_values = obj.retrieveParameterValues(combos{p,:}, 1:3);
                new_rows = [files, repmat(tab_values, n_files, 1)];
                tab{end+1:end+n_files, :} = new_rows;
            end

            if height(tab) == 0
                warning("DatasetFileManager:noFilesFound", ...
                    "No files found matching pattern '%s'. Check your OriginalFilepathPattern and parameter labels.", ...
                    obj.OriginalFilepathPattern);
            end
        end

        function disp(obj)
            %Display DATASETFILEMANAGER instance
            fprintf("DatasetFileManager object\n\n")
            fprintf("\t OriginalFilepathPattern: ""%s""\n", obj.OriginalFilepathPattern)
            fprintf("\tExtractedFilepathPattern: ""%s""\n", obj.ExtractedFilepathPattern)
            for p = 1:length(obj.ParameterTables)
                fprintf("\nParameter ""%s"":\n\n", obj.ParameterNames(p));
                disp(obj.ParameterTables{p})
            end
        end
    end

    methods (Access=private)
        function ind = findParameter(obj, par_name)
            %Retrieve parameter index
            ind = find(obj.ParameterNames == par_name, 1);
        end

        function [par_ind, label_ind] = findlabel(obj, par_name, label)
            %Retrieve parameter expression index from parameter table
            par_ind = obj.findParameter(par_name);
            if ~isempty(par_ind)
                label_ind = find(obj.ParameterTables{par_ind}{:,1} == label);
            else
                label_ind = [];
            end
        end

        function combos = labelIndexCombos(obj)
            %Create all parameter expression index combinations
            table_row_inds = @(tab) 1:height(tab);
            label_inds = cellfun(table_row_inds, obj.ParameterTables, UniformOutput=false);
            combos = combinations(label_inds{:});
        end

        function tab = initializeFileTable(obj)
            %Create empty file table (used in compileOriginalFileTable())
            n_pars = length(obj.ParameterTables);
            t_size = [0, 3*n_pars+1];
            t_types = repmat("string", 1, t_size(2));
            table_var_names = @(tab) string(tab.Properties.VariableNames);
            t_names = ["Filepath", cell2mat(cellfun(table_var_names, obj.ParameterTables', UniformOutput=false))];
            tab = table( ...
                Size=t_size, ...
                VariableTypes=t_types, ...
                VariableNames=t_names);
        end

        function vals = retrieveParameterValues(obj, row_selectors, cols)
            %Retrieve combination of specific parameter expressions
            n_pars = length(obj.ParameterTables);
            n_cols = length(cols);
            assert(length(row_selectors) == n_pars)
            vals = strings(1, n_pars * n_cols);
            for p = 1:n_pars
                ind1 = (p-1) * n_cols + 1;
                ind2 = p * n_cols;
                vals(ind1:ind2) = obj.ParameterTables{p}{row_selectors(p), cols};
            end
        end
    end

    methods (Static, Access=private)
        function tab = initializeParameterTable(name)
            %Create empty parameter table (used in addParameter())
            vars = [ ...
                compose("%sLabel", name)                     , "string"; ...
                compose("%sOriginalFilepathComponent", name) , "string"; ...
                compose("%sExtractedFilepathComponent", name), "string"; ...
                ];
            tabSize = [0, 3];
            tab = table(Size=tabSize, VariableNames=vars(:,1), VariableTypes=vars(:,2));
        end
    end
end