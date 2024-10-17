% This MATLAB code depicts the construction of a phylogenetic tree.
% It allows building the phylogenetic tree using the Neighbor-Joining 
% and UPGMA methods.

% The person executing this code first needs to choose what type of 
% sequences will be used for this analysis, whether they are sequences to be 
% aligned or already aligned sequences.

% If they are sequences to be aligned, it is necessary to provide a text 
% document containing the names and GenBank IDs in this format:
% [name_corresponding_to_sequence     GenBank ID]
% ...
% The sequences will then be downloaded and aligned according to a 
% multi-align and the scoringMatrix is nuc44.

% If they are already aligned, they must be in a fasta file 
% that will be read from which the name and sequence will be extracted.

% In the next step, the user will be able to select what type of 
% sequence they will analyze. The process is the same for both options; the only 
% difference is the preparation of the sequences. 
% The code is divided this way to avoid problems.

disp('some steps may take time')
%% BEGINNING OF THE CODE

% Main menu for option selection
disp('Option 1 ----> sequences to be aligned');
disp('Option 2 ----> already aligned sequences');
option = input('Select an option (1 or 2): ');

if option == 1
    %% PREPARATION OF SEQUENCES
    % Write the name of the .txt file that contains the information
    ficheiro = input('Name of the file with the information: ', 's'); % Requests the file name

    type(ficheiro); % Displays the content of the file

    dados = readcell(ficheiro); % Reads the data from the file

    % Passes the GenBank IDs of the sequences to the variable 'ids'
    ids = dados(:, 2); % Column containing the GenBank IDs

    % The sequences will be downloaded

    num_seqs = length(ids); % Creates a variable that indicates the number of sequences
    sequencias = cell(num_seqs, 1); % Creates a variable with the required number of cells

    % Download
    for i = 1:num_seqs
        accession = ids{i};
        seqData = getgenbank(accession);
        sequencias{i} = seqData.Sequence;
    end

    % Displays the names of the downloaded sequences
    nomes = dados(:, 1); % Column containing the names

    %% Alignment of sequences

    % Converts the downloaded sequences into a format suitable for multiple alignment
    seqArray = cellfun(@(x) x, sequencias, 'UniformOutput', false);
    seqArray = seqArray(~cellfun('isempty', seqArray));

    % Perform multiple alignment with the ScoringMatrix nuc44
    alinhamento = multialign(seqArray, 'ScoringMatrix', 'nuc44');

    % Creates a variable to store the aligned sequences
    sequencias_alinhadas = cell(num_seqs, 1);

    % Stores the aligned sequences 
    for i = 1:num_seqs
        sequencias_alinhadas{i} = alinhamento(i, :);
    end

    %% CONSTRUCTION OF TREES
    % Here you will choose what type of tree you want

    % Main menu for selecting the type of tree to create 
    disp('Option 1 ----> UPGMA phylogenetic tree');
    disp('Option 2 ----> Neighbor-Joining phylogenetic tree');
    option = input('Select an option (1 or 2): ');

    % Creation of an if/else code, which will be used 
    % to execute the corresponding code for each tree
    if option == 1
        % UPGMA Tree
        disp('UPGMA tree');
        % Creates a UPGMA tree, the distance is calculated using the Jukes-Cantor method

        % The evolutionary distance is calculated
        distances = seqpdist(sequencias_alinhadas, 'Method', 'Jukes-Cantor'); 

        % Builds the UPGMA tree using the calculated distances
        UPGMAtree = seqlinkage(distances, 'UPGMA', nomes);

        % Bootstrap for evaluating the tree's confidence
        num_boots = 100;
        seq_len = length(sequencias_alinhadas{1});
        
        % Initializes cells to store bootstrap sequences
        boots = cell(num_boots, 1);

        % Generates bootstrap replicates of the sequences
        for n = 1:num_boots
            reorder_index = randsample(seq_len, seq_len, true); % Generates random indices
            bootseq = struct('Header', nomes, 'Sequence', sequencias_alinhadas);
            % Reorders the sequences according to the random indices
            for i = 1:num_seqs
                % Reorders sequence i based on the random indices
                bootseq(i).Sequence = sequencias_alinhadas{i}(reorder_index);
            end
            boots{n} = bootseq; % Stores the bootstrap replicate
        end

        % Anonymous function to construct the tree using UPGMA and calculate the distance
        fun = @(x) seqlinkage(x, 'average', nomes);
        boot_trees = cell(num_boots, 1); % Initializes cells to store bootstrap trees

        % Calculates the distance between bootstrap replicates and builds the corresponding trees
        for n = 1:num_boots
            dist_tmp = seqpdist({boots{n}.Sequence});
            boot_trees{n} = fun(dist_tmp);
        end

        % Evaluation of the topology of the original tree concerning the bootstrap replicates
        orig_pointers = cell(num_seqs-1, 1);
        orig_species = cell(num_seqs-1, 1);

        % Obtains the topology and species for each branch of the original tree
        for i = num_seqs-1:-1:1
            branch_pointer = i + num_seqs;
            sub_tree = subtree(UPGMAtree, branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree, 'LeafNames'));
        end

        % Evaluation of the topology of the bootstrap trees
        clusters_pointers = cell(num_seqs-1, num_boots);
        clusters_species = cell(num_seqs-1, num_boots);

        % Obtains the topology and species for each branch of the bootstrap trees
        for j = 1:num_boots
            for i = num_seqs-1:-1:1
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j}, branch_ptr);
                clusters_pointers{i, j} = getcanonical(sub_tree);
                clusters_species{i, j} = sort(get(sub_tree, 'LeafNames'));
            end
        end
        
        % Counting branches with similar topology between the original tree and bootstrap trees
        count = zeros(num_seqs-1, 1);
        for i = 1:num_seqs-1
            for j = 1:num_boots
                if isequal(orig_pointers{i}, clusters_pointers{i, j}) && ...
                   isequal(orig_species{i}, clusters_species{i, j})
                    count(i) = count(i) + 1;
                end
            end
        end
        
        % Calculation of confidence probability (Pc)
        Pc = count ./ num_boots;   % confidence (Pc)

        % Visualizing the Confidence Values in the Original Tree
        [ptrs, dist, names] = get(UPGMAtree, 'POINTERS', 'DISTANCES', 'NODENAMES');
        for i = 1:num_seqs-1
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr}, ', confidence: ', num2str(100*Pc(i)), ' %'];
        end

        % Builds the final tree
        tr = phytree(ptrs, dist, names);

        % Shows the branch points with a confidence level greater than C
        % A higher confidence level (e.g., 99%) may provide a more accurate estimate, 
        % A lower confidence level (e.g., 90%) may result 
        % in wider and less precise confidence intervals.
        c = input('Confidence level greater than: ');
        if c >= 0 && c <= 1
            high_conf_branch_ptr = find(Pc > c) + num_seqs;
            view(tr, high_conf_branch_ptr);
        else
            disp('Please enter the confidence level again, invalid value.');
        end

        disp('Process completed successfully.');
        
    elseif option == 2
        % Neighbor-Joining Tree
        disp('Neighbor-Joining tree');
        % Calculates the distances between the aligned sequences using the Jukes-Cantor method
        distances = seqpdist(sequencias_alinhadas, 'Method', 'Jukes-Cantor');
        % Builds the Neighbor-Joining tree using the calculated distances
        NJtree = seqneighjoin(distances, 'equivar', nomes);

        % Bootstrap for evaluating the tree's confidence
        num_boots = 100;
        seq_len = length(sequencias_alinhadas{1});
        
        boots = cell(num_boots, 1); % Initializes cells to store bootstrap sequences

        % Generates bootstrap replicates of the sequences
        for n = 1:num_boots
            reorder_index = randsample(seq_len, seq_len, true); % Generates random indices
            bootseq = struct('Header', nomes, 'Sequence', sequencias_alinhadas);
            % Reorders the sequences according to the random indices
            for i = 1:num_seqs
                bootseq(i).Sequence = sequencias_alinhadas{i}(reorder_index);
            end
            boots{n} = bootseq; % Stores the bootstrap replicate
        end

        % Anonymous function to construct the tree using Neighbor-Joining and calculate the distance
        fun = @(x) seqneighjoin(x, 'equivar', nomes);
        boot_trees = cell(num_boots, 1); % Initializes cells to store bootstrap trees

        % Calculates the distance between bootstrap replicates and builds the corresponding trees
        for n = 1:num_boots
            dist_tmp = seqpdist({boots{n}.Sequence});
            boot_trees{n} = fun(dist_tmp);
        end

        % Evaluation of the topology of the original tree concerning the bootstrap replicates
        orig_pointers = cell(num_seqs-1, 1);
        orig_species = cell(num_seqs-1, 1);

        % Obtains the topology and species for each branch of the original tree
        for i = num_seqs-1:-1:1
            branch_pointer = i + num_seqs;
            sub_tree = subtree(NJtree, branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree, 'LeafNames'));
        end

        % Evaluation of the topology of the bootstrap trees
        clusters_pointers = cell(num_seqs-1, num_boots);
        clusters_species = cell(num_seqs-1, num_boots);

        % Obtains the topology and species for each branch of the bootstrap trees
        for j = 1:num_boots
            for i = num_seqs-1:-1:1
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j}, branch_ptr);
                clusters_pointers{i, j} = getcanonical(sub_tree);
                clusters_species{i, j} = sort(get(sub_tree, 'LeafNames'));
            end
        end
        
        % Counting branches with similar topology between the original tree and bootstrap trees
        count = zeros(num_seqs-1, 1);
        for i = 1:num_seqs-1
            for j = 1:num_boots
                if isequal(orig_pointers{i}, clusters_pointers{i, j}) && ...
                   isequal(orig_species{i}, clusters_species{i, j})
                    count(i) = count(i) + 1;
                end
            end
        end
        
        % Calculation of confidence probability (Pc)
        Pc = count ./ num_boots;   % confidence (Pc)

        % Visualizing the Confidence Values in the Original Tree
        [ptrs, dist, names] = get(NJtree, 'POINTERS', 'DISTANCES', 'NODENAMES');
        for i = 1:num_seqs-1
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr}, ', confidence: ', num2str(100*Pc(i)), ' %'];
        end

        % Builds the final tree
        tr = phytree(ptrs, dist, names);

        % Shows the branch points with a confidence level greater than C
        % A higher confidence level (e.g., 99%) may provide a more accurate estimate, 
        % A lower confidence level (e.g., 90%) may result 
        % in wider and less precise confidence intervals.
        c = input('Confidence level greater than: ');
        if c >= 0 && c <= 1
            high_conf_branch_ptr = find(Pc > c) + num_seqs;
            view(tr, high_conf_branch_ptr);
        else
            disp('Please enter the confidence level again, invalid value.');
        end

        disp('Process completed successfully.');
        
    end

elseif option == 2
    %% READING THE ALIGNED SEQUENCES
    
    primates = fastaread('primatesaligned.fa');
    num_seqs = length(primates);
    % Exibe o nome de cada sequencia
    for i = 1:num_seqs
        fprintf('Nome_da_Sequencia: %s\n', primates(i).Header);     
    end
    
    %% TREE CONSTRUCTION

    % Main menu for selecting the type of tree created 
    % Ask the user which option they want to select
    disp('Option 1 ----> UPGMA phylogenetic tree');
    disp('Opção 2 ----> árvore filogenética Neighbor-Joining');
    option = input('Select an option (1 or 2):  ');
    
    
    % Creation of an if/else code, which will be used to execute
    % the code corresponding to each tree
    if option == 1
        disp('UPGMA tree');
        % Calculates the distances between sequences using the Jukes-Cantor method
        distances = seqpdist(primates, 'Method', 'Jukes-Cantor');
        % Constructs the UPGMA phylogenetic tree using the calculated distances
        UPGMAtree = seqlinkage(distances, 'UPGMA', primates);
    
    %% BOOTSTRAPING----------------------------------------------------------
              
        % Bootstrap for tree confidence assessment
        num_boots = 100; % Number of bootstrap repetitions
        seq_len = length(primates(1).Sequence); % Comprimento das sequências
        
        % Initializes cells to store bootstrap sequences
        boots = cell(num_boots,1);

        % Generates bootstrap replicates of the sequences
        for n = 1:num_boots
            reorder_index = randsample(seq_len,seq_len,true); % Generates random indexes
            for i = num_seqs:-1:1 %reverse order for pre-allocation  
                % Copies the original sequence header
                bootseq(i).Header = primates(i).Header;
                % Reorders the sequence based on random indexes and removes gaps
                bootseq(i).Sequence = strrep(primates(i).Sequence(reorder_index),'-',''); 
            end
            % Stores the bootstrap replicate
            boots{n} = bootseq; 
        end
        
        % Anonymous function to build the tree using UPGMA and calculate the distance
        fun = @(x) seqlinkage(x,'average',{primates.Header});

        % Initializes cells to store bootstrap trees
        boot_trees = cell(num_boots,1); 
       
          
        % Calculates the distance between bootstrap replicates and builds the corresponding trees
        for n = 1:num_boots
            dist_tmp = seqpdist(boots{n});
            boot_trees{n} = fun(dist_tmp);
        end
    
    
        % Obtains the topology and species for each branch of the original tree
        for i = num_seqs-1:-1:1  % for every branch, reverse order to preallocate
            branch_pointer = i + num_seqs;
            sub_tree = subtree(UPGMAtree,branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree,'LeafNames'));
        end


        % Obtains the topology and species for each branch of the bootstrap trees
        for j = num_boots:-1:1
            for i = num_seqs-1:-1:1  % for every branch
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j},branch_ptr);
                clusters_pointers{i,j} = getcanonical(sub_tree);
                clusters_species{i,j} = sort(get(sub_tree,'LeafNames'));
            end
        end

        % Count of branches with similar topology between the original tree and the bootstrap trees
        count = zeros(num_seqs-1,1);
        for i = 1 : num_seqs -1  % for every branch
            for j = 1 : num_boots * (num_seqs-1)
                if isequal(orig_pointers{i},clusters_pointers{j})
                    if isequal(orig_species{i},clusters_species{j})
                        count(i) = count(i) + 1;
                    end
                end
            end
        end
        
        % Calculation of the probability of confidence (Pc)
        Pc = count ./ num_boots;   % confidence probability (Pc)
    
    
        %Visualizing the Confidence Values in the Original Tree
        [ptrs,dist,names] = get(UPGMAtree,'POINTERS','DISTANCES','NODENAMES');
        
        for i = 1:num_seqs -1  % for every branch
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
        end
    
     
        % Building the final tree   
        tr = phytree(ptrs,dist,names); 
    

        % Displays the bifurcation branches with a confidence level higher than X
        % A higher confidence level (e.g. 99%) may provide a more accurate estimate, 
        % A lower confidence level (e.g. 90%) may result in 
        % in wider and less precise confidence intervals.
        
        c = input('confidence level greater than: ');
        
        % Checking if c is within the range
        if c >= 0 && c <= 1
             high_conf_branch_ptr = find(Pc > c) + num_seqs;
             view(tr, high_conf_branch_ptr);
        else
            disp('Re-enter the confidence level, invalid value.');
        end

        disp('Process completed successfully.');
    
    elseif option == 2
        disp('Neighbor-Joining tree');
        % Calculates the distances between aligned sequences using the Jukes-Cantor method
        distances = seqpdist(primates, 'Method', 'Jukes-Cantor');
        % Builds the Neighbor-Joining tree using the calculated distances
        NJtree = seqneighjoin(distances, 'equivar', primates);
    
        %% BOOTSTRAPING----------------------------------------------------------
        
        % Bootstrap for tree confidence assessment    
        num_boots = 100; % Number of bootstrap repetitions
        seq_len = length(primates(1).Sequence); % Length of sequences
          
        % Initializes cells to store bootstrap sequences
        boots = cell(num_boots,1);

        % Generates bootstrap replicates of the sequences
        for n = 1:num_boots
            reorder_index = randsample(seq_len,seq_len,true); % Generates random indexes
            for i = num_seqs:-1:1 %reverse order to preallocate memory
                % Copies the original sequence header
                bootseq(i).Header = primates(i).Header;
                % Reorders the sequence based on random indexes and removes gaps
                bootseq(i).Sequence = strrep(primates(i).Sequence(reorder_index),'-',''); 
            end

            % Stores the bootstrap replicate
            boots{n} = bootseq; 
        end
        
        % Anonymous function to build the tree using Neighbor-Joining and calculate the distance
        fun = @(x) seqlinkage(x,'average',{primates.Header});

        % Initializes cells to store bootstrap trees
        boot_trees = cell(num_boots,1); 
       
             
        % Calculates the distance between bootstrap replicates and builds the corresponding trees
        for n = 1:num_boots            
            dist_tmp = seqpdist(boots{n});
            boot_trees{n} = fun(dist_tmp);
        end
    
        
        % Obtains the topology and species for each branch of the original tree
        for i = num_seqs-1:-1:1  % for every branch, reverse order to preallocate
            branch_pointer = i + num_seqs;
            sub_tree = subtree(NJtree,branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree,'LeafNames'));
        end
        
        % Obtains the topology and species for each branch of the bootstrap trees
        for j = num_boots:-1:1
            for i = num_seqs-1:-1:1  % for every branch
                branch_ptr = i + num_seqs;
                % should be this line
                sub_tree = subtree(boot_trees{j},branch_ptr);
                clusters_pointers{i,j} = getcanonical(sub_tree);
                clusters_species{i,j} = sort(get(sub_tree,'LeafNames'));
            end
        end
       
        % Count of branches with similar topology between the original tree and the bootstrap trees
        count = zeros(num_seqs-1,1);
        for i = 1 : num_seqs -1  % for every branch
            for j = 1 : num_boots * (num_seqs-1)
                if isequal(orig_pointers{i},clusters_pointers{j})
                    if isequal(orig_species{i},clusters_species{j})
                        count(i) = count(i) + 1;
                    end
                end
            end
        end

        % Calculation of the probability of confidence (Pc)
        Pc = count ./ num_boots;   % confidence probability (Pc)
    
        % Display of confidence values in the original tree
        [ptrs,dist,names] = get(NJtree,'POINTERS','DISTANCES','NODENAMES');
        
        for i = 1:num_seqs -1  % for every branch
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
        end
    
        % Build the final tree   
        tr = phytree(ptrs,dist,names); 
    

        % Show the branch points that have a confidence level greater than C
        % A higher confidence level (e.g. 99%) may provide a more accurate estimate, 
        % A lower confidence level (e.g. 90%) may result 
        % in wider and less precise confidence intervals.          
        c = input('confidence level greater than: ');
        
        % Checking if c is within the range
        if c >= 0 && c <= 1
             high_conf_branch_ptr = find(Pc > c) + num_seqs;
             view(tr, high_conf_branch_ptr);
        else
            disp('Re-enter the confidence level, invalid value.');
        end

        disp('Processo concluído com sucesso.');
    end
end

disp('')
