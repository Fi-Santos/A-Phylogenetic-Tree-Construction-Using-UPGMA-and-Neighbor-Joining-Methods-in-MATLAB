% Este codigo de MATLAB retrata a construção de uma arvore filogenetica.
% Ele premite construir a arvore filogenetica utilizando os metodos Neighbor-joining
% e UPGMA

% A pessoa que executa este codigo primeiro tem de escolher que tipo de
% sequencias vão ser utilizadas para esta analize, se são sequencias por alinhar
% ou se são sequencias ja alinhadas.

% Se forem sequencias por alinhar é preciso fornecer um documento texto onde
% contem os nomes e o GenBank IDs , neste formato:
% [nome_correspondente_a_sequencia     GenBank ID]
% ...
% É então feito o download das sequencias e são alinhadas segundo um
% multialign e a scoringMatrix é a nuc44

% Se elas já estiverem alinhadas , elas têm de estar num ficheiro fasta 
% que vai ser lido e de onde se vai retirar o nome e a sequencia 


% Na sequeção seguinte o utilizador vai poder selecionar que tipo de sequencia vai analizar
% O processo é o mesmo para as duas opções a unica coisa diferente é a
% preparação das sequencias 
% O codigo esta dividido desta forma para evitar problemas

disp('algumas etapas podem demorar')
%% COMEÇO DO CODIGO
 
% Menu principal para seleção de opções
disp('Opção 1 ----> sequencias por alinhar');
disp('Opção 2 ----> sequencias já alinhadas');
option = input('Selecione uma opção (1 ou 2): ');

if option == 1
    %% PREPARAÇÃO DAS SEQUÊNCIAS
    % Escrever o  nome do ficheiro.txt que contem as informações
    ficheiro = input('Nome do ficheiro com as informações: ', 's'); % Solicita o nome do ficheiro

    type(ficheiro); % Exibe o conteúdo do ficheiro

    dados = readcell(ficheiro); % Le os dados do ficheiro

    %Passa para a variavel 'ids' os GenBankIDs das sequencias
    ids = dados(:, 2); % Coluna contendo os IDs do GenBank

    % Vai ser feito o Download das sequências

    num_seqs = length(ids); % Cria uma variavel que diz o numero de sequencias
    sequencias = cell(num_seqs, 1); % Cria uma varivel com o numero de celulas necessario 
    
    % Download
    for i = 1:num_seqs
        accession = ids{i};
        seqData = getgenbank(accession);
        sequencias{i} = seqData.Sequence;
    end

    % Exibe o nome das sequências baixadas
    nomes = dados(:, 1); % Coluna contendo os nomes

    %% Alinhamento das sequencias
    
    % Converte as sequências baixadas em um formato adequado para multialinhamento
    seqArray = cellfun(@(x) x, sequencias, 'UniformOutput', false);
    seqArray = seqArray(~cellfun('isempty', seqArray));

    % Realizar o alinhamento múltiplo com a ScoringMatrix nuc44
    alinhamento = multialign(seqArray, 'ScoringMatrix', 'nuc44');

    % Cria uma variável para armazenar as sequências alinhadas
    sequencias_alinhadas = cell(num_seqs, 1);
    
    % Armazena as sequencias alinhadas 
    for i = 1:num_seqs
        sequencias_alinhadas{i} = alinhamento(i, :);
    end


    %% CONSTRUÇÃO DAS ÁRVORES
    % Aqui vai-se escolher que tipo de arvore se quer

    % Menu principal para seleção do tipo de arvore criada 
    disp('Opção 1 ----> árvore filogenética UPGMA');
    disp('Opção 2 ----> árvore filogenética Neighbor-Joining');
    option = input('Selecione uma opção (1 ou 2): ');

    % Criação de um codigo if/else , vai ser usadon para que seja executado
    % o codigo correspondente a cada arvore
    if option == 1
        % Árvore UPGMA
        disp('UPGMA tree');
        % Cria uma arvore UPGMA , a distancia é calculada usando o metodo Jukes-Cantor

        % É calculada a distancia evolutiva
        distances = seqpdist(sequencias_alinhadas, 'Method', 'Jukes-Cantor'); 

        % Constrói a árvore UPGMA usando as distâncias calculadas
        UPGMAtree = seqlinkage(distances, 'UPGMA', nomes);

        % Bootstrap para avaliação da confiança da árvore
        num_boots = 100;
        seq_len = length(sequencias_alinhadas{1});
        
        % Inicializa células para armazenar as sequências bootstrap
        boots = cell(num_boots, 1);

        % Gera replicatas bootstrap das sequências
        for n = 1:num_boots
            reorder_index = randsample(seq_len, seq_len, true); % Gera índices aleatórios
            bootseq = struct('Header', nomes, 'Sequence', sequencias_alinhadas);
            % Reordena as sequências de acordo com os índices aleatórios
            for i = 1:num_seqs
                % Reordena a sequência i com base nos índices aleatórios
                bootseq(i).Sequence = sequencias_alinhadas{i}(reorder_index);
            end
            boots{n} = bootseq; % Armazena a replicata bootstrap
        end

        %Função anônima para construir a árvore usando UPGMA e calcular a distância
        fun = @(x) seqlinkage(x, 'average', nomes);
        boot_trees = cell(num_boots, 1); % Inicializa células para armazenar as árvores bootstrap

        % Calcula a distância entre as replicatas bootstrap e constrói as árvores correspondentes
        for n = 1:num_boots
            dist_tmp = seqpdist({boots{n}.Sequence});
            boot_trees{n} = fun(dist_tmp);
        end

        % Avaliação da topologia da árvore original em relação às replicatas bootstrap
        orig_pointers = cell(num_seqs-1, 1);
        orig_species = cell(num_seqs-1, 1);

        % Obtém a topologia e as espécies para cada ramo da árvore original
        for i = num_seqs-1:-1:1
            branch_pointer = i + num_seqs;
            sub_tree = subtree(UPGMAtree, branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree, 'LeafNames'));
        end

        % Avaliação da topologia das árvores bootstrap
        clusters_pointers = cell(num_seqs-1, num_boots);
        clusters_species = cell(num_seqs-1, num_boots);

        % Obtém a topologia e as espécies para cada ramo das árvores bootstrap
        for j = 1:num_boots
            for i = num_seqs-1:-1:1
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j}, branch_ptr);
                clusters_pointers{i, j} = getcanonical(sub_tree);
                clusters_species{i, j} = sort(get(sub_tree, 'LeafNames'));
            end
        end
        
        % Contagem de ramos com topologia similar entre a árvore original e as árvores bootstrap
        count = zeros(num_seqs-1, 1);
        for i = 1:num_seqs-1
            for j = 1:num_boots
                if isequal(orig_pointers{i}, clusters_pointers{i, j}) && ...
                   isequal(orig_species{i}, clusters_species{i, j})
                    count(i) = count(i) + 1;
                end
            end
        end
        
        % Cálculo da probabilidade de confiança (Pc)
        Pc = count ./ num_boots;   % confiança (Pc)

        % Visualizar os Valores de Confiança na Árvore Original
        [ptrs, dist, names] = get(UPGMAtree, 'POINTERS', 'DISTANCES', 'NODENAMES');
        for i = 1:num_seqs-1
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr}, ', confiança: ', num2str(100*Pc(i)), ' %'];
        end

        % Constrói a árvore final
        tr = phytree(ptrs, dist, names);

        % Mostra os pontos de ramificação com um nível de confiança maior que C
        % Um nível de confiança mais alto (por exemplo, 99%) possa fornecer uma estimativa mais precisa, 
        % Um nível de confiança mais baixo (por exemplo, 90%) pode resultar 
        % em intervalos de confiança mais amplos e menos precisos.
        c = input('Nível de confiança maior que: ');
        if c >= 0 && c <= 1
            high_conf_branch_ptr = find(Pc > c) + num_seqs;
            view(tr, high_conf_branch_ptr);
        else
            disp('Insira novamente o nível de confiança, valor inválido.');
        end

        disp('Processo concluído com sucesso.');
        
    elseif option == 2
        % Árvore Neighbor-Joining
        disp('Neighbor-Joining tree');
        % Calcula as distâncias entre as sequências alinhadas usando o método de Jukes-Cantor
        distances = seqpdist(sequencias_alinhadas, 'Method', 'Jukes-Cantor');
        % Constrói a árvore Neighbor-Joining usando as distâncias calculadas
        NJtree = seqneighjoin(distances, 'equivar', nomes);

        % Bootstrap para avaliação da confiança da árvore
        num_boots = 100;
        seq_len = length(sequencias_alinhadas{1});
        
        boots = cell(num_boots, 1); % Inicializa células para armazenar as sequências bootstrap

        % Gera replicatas bootstrap das sequências
        for n = 1:num_boots
            reorder_index = randsample(seq_len, seq_len, true); % Gera índices aleatórios
            bootseq = struct('Header', nomes, 'Sequence', sequencias_alinhadas);

            % Reordena as sequências de acordo com os índices aleatórios
            for i = 1:num_seqs
                % Reordena a sequência i com base nos índices aleatórios
                bootseq(i).Sequence = sequencias_alinhadas{i}(reorder_index);
            end
            % Armazena a replicata bootstrap
            boots{n} = bootseq;
        end

        % Função anônima para construir a árvore usando Neighbor-Joining e calcular a distância
        fun = @(x) seqlinkage(x, 'average', nomes);
        % Inicializa células para armazenar as árvores bootstrap
        boot_trees = cell(num_boots, 1);

        % Calcula a distância entre as replicatas bootstrap e constrói as árvores correspondentes
        for n = 1:num_boots
            dist_tmp = seqpdist({boots{n}.Sequence});
            boot_trees{n} = fun(dist_tmp);
        end

        % Avaliação da topologia da árvore original em relação às replicatas bootstrap
        orig_pointers = cell(num_seqs-1, 1);
        orig_species = cell(num_seqs-1, 1);

        % Obtém a topologia e as espécies para cada ramo da árvore original
        for i = num_seqs-1:-1:1
            branch_pointer = i + num_seqs;
            sub_tree = subtree(NJtree, branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree, 'LeafNames'));
        end
        
        % Avaliação da topologia das árvores bootstrap
        clusters_pointers = cell(num_seqs-1, num_boots);
        clusters_species = cell(num_seqs-1, num_boots);

        % Obtém a topologia e as espécies para cada ramo das árvores bootstrap
        for j = 1:num_boots
            for i = num_seqs-1:-1:1
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j}, branch_ptr);
                clusters_pointers{i, j} = getcanonical(sub_tree);
                clusters_species{i, j} = sort(get(sub_tree, 'LeafNames'));
            end
        end
        
        % Contagem de ramos com topologia similar entre a árvore original e as árvores bootstrap
        count = zeros(num_seqs-1, 1);
        for i = 1:num_seqs-1
            for j = 1:num_boots
                if isequal(orig_pointers{i}, clusters_pointers{i, j}) && ...
                   isequal(orig_species{i}, clusters_species{i, j})
                    count(i) = count(i) + 1;
                end
            end
        end
        

        % Cálculo da probabilidade de confiança (Pc)
        Pc = count ./ num_boots;   % confiança (Pc)

        % Visualizando os Valores de Confiança na Árvore Original
        [ptrs, dist, names] = get(NJtree, 'POINTERS', 'DISTANCES', 'NODENAMES');
        for i = 1:num_seqs-1
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr}, ', confiança: ', num2str(100*Pc(i)), ' %'];
        end

        % Constrói a árvore final
        tr = phytree(ptrs, dist, names);

        % Mostrar os pontos de ramificação com um nível de confiança maior que C
        % Um nível de confiança mais alto (por exemplo, 99%) possa fornecer uma estimativa mais precisa, 
        % Um nível de confiança mais baixo (por exemplo, 90%) pode resultar 
        % em intervalos de confiança mais amplos e menos precisos.

        c = input('Nível de confiança maior que: ');
        if c >= 0 && c <= 1
            high_conf_branch_ptr = find(Pc > c) + num_seqs;
            view(tr, high_conf_branch_ptr);
        else
            disp('Insira novamente o nível de confiança, valor inválido.');
        end

        disp('Processo concluído com sucesso.');

    else
        disp('Opção inválida, por favor selecione 1 ou 2.');
    end

elseif option == 2
    %% LEITURA DAS SEQUÊNCIAS ALINHADAS
    
    primates = fastaread('primatesaligned.fa');
    num_seqs = length(primates);
    % Exibe o nome de cada sequencia
    for i = 1:num_seqs
        fprintf('Nome_da_Sequencia: %s\n', primates(i).Header);     
    end

    %% CONSTRUÇÃO DAS ARVORES

    % Menu principal para seleção do tipo de arvore criada 
    % Perguntar ao usuário qual opção deseja selecionar
    disp('Opção 1 ----> árvore filogenética UPGMA');
    disp('Opção 2 ----> árvore filogenética Neighbor-Joining');
    option = input('Selecione uma opção (1 ou 2): ');
    
    
    % Criação de um codigo if/else , vai ser usadon para que seja executado
    % o codigo correspondente a cada arvore
    if option == 1
        disp('UPGMA tree');
        % Calcula as distâncias entre as sequências usando o método de Jukes-Cantor
        distances = seqpdist(primates, 'Method', 'Jukes-Cantor');
        % Constrói a árvore filogenética UPGMA usando as distâncias calculadas
        UPGMAtree = seqlinkage(distances, 'UPGMA', primates);
    
    %% BOOTSTRAPING----------------------------------------------------------
              
        % Bootstrap para avaliação da confiança da árvore
        num_boots = 100; % Número de repetições de bootstrap
        seq_len = length(primates(1).Sequence); % Comprimento das sequências
        
        % Inicializa células para armazenar as sequências bootstrap
        boots = cell(num_boots,1);

        % Gera replicatas bootstrap das sequências
        for n = 1:num_boots
            reorder_index = randsample(seq_len,seq_len,true); % Gera índices aleatórios
            for i = num_seqs:-1:1 %reverse order to preallocate 
                % Copia o cabeçalho da sequência original
                bootseq(i).Header = primates(i).Header;
                % Reordena a sequência com base nos índices aleatórios e remove os gaps
                bootseq(i).Sequence = strrep(primates(i).Sequence(reorder_index),'-',''); 
            end
            % Armazena a replicata bootstrap
            boots{n} = bootseq; 
        end
        
        % Função anônima para construir a árvore usando UPGMA e calcular a distância
        fun = @(x) seqlinkage(x,'average',{primates.Header});

        % Inicializa células para armazenar as árvores bootstrap
        boot_trees = cell(num_boots,1); 
       
          
        % Calcula a distância entre as replicatas bootstrap e constrói as árvores correspondentes
        for n = 1:num_boots
            dist_tmp = seqpdist(boots{n});
            boot_trees{n} = fun(dist_tmp);
        end
    
    
        % Obtém a topologia e as espécies para cada ramo da árvore original
        for i = num_seqs-1:-1:1  % for every branch, reverse order to preallocate
            branch_pointer = i + num_seqs;
            sub_tree = subtree(UPGMAtree,branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree,'LeafNames'));
        end


        % Obtém a topologia e as espécies para cada ramo das árvores bootstrap
        for j = num_boots:-1:1
            for i = num_seqs-1:-1:1  % for every branch
                branch_ptr = i + num_seqs;
                sub_tree = subtree(boot_trees{j},branch_ptr);
                clusters_pointers{i,j} = getcanonical(sub_tree);
                clusters_species{i,j} = sort(get(sub_tree,'LeafNames'));
            end
        end

        % Contagem de ramos com topologia similar entre a árvore original e as árvores bootstrap
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
        
        % Cálculo da probabilidade de confiança (Pc)
        Pc = count ./ num_boots;   % confidence probability (Pc)
    
    
        %Visualizing the Confidence Values in the Original Tree
        [ptrs,dist,names] = get(UPGMAtree,'POINTERS','DISTANCES','NODENAMES');
        
        for i = 1:num_seqs -1  % for every branch
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
        end
    
     
        % Construir a árvore final   
        tr = phytree(ptrs,dist,names); 
    

        % Exibe os ramos de bifurcação com um nível de confiança maior que X
        % Um nível de confiança mais alto (por exemplo, 99%) possa fornecer uma estimativa mais precisa, 
        % Um nível de confiança mais baixo (por exemplo, 90%) pode resultar 
        % em intervalos de confiança mais amplos e menos precisos.
        
        c = input('nivel de confiança maior que: ');
        
        % Verificando se c está dentro do intervalo
        if c >= 0 && c <= 1
             high_conf_branch_ptr = find(Pc > c) + num_seqs;
             view(tr, high_conf_branch_ptr);
        else
            disp('Insira novamente o nível de confiança, valor inválido.');
        end

        disp('Processo concluído com sucesso.');
    
    elseif option == 2
        disp('Neighbor-Joining tree');
        % Calcula as distâncias entre as sequências alinhadas usando o método de Jukes-Cantor
        distances = seqpdist(primates, 'Method', 'Jukes-Cantor');
        % Constrói a árvore Neighbor-Joining usando as distâncias calculadas
        NJtree = seqneighjoin(distances, 'equivar', primates);
    
        %% BOOTSTRAPING----------------------------------------------------------
        
        % Bootstrap para avaliação da confiança da árvore    
        num_boots = 100; % Número de repetições de bootstrap
        seq_len = length(primates(1).Sequence); % Comprimento das sequências
          
        % Inicializa células para armazenar as sequências bootstrap
        boots = cell(num_boots,1);

        % Gera replicatas bootstrap das sequências
        for n = 1:num_boots
            reorder_index = randsample(seq_len,seq_len,true); % Gera índices aleatórios
            for i = num_seqs:-1:1 %reverse order to preallocate memory
                % Copia o cabeçalho da sequência original
                bootseq(i).Header = primates(i).Header;
                % Reordena a sequência com base nos índices aleatórios e remove os gaps
                bootseq(i).Sequence = strrep(primates(i).Sequence(reorder_index),'-',''); 
            end

            % Armazena a replicata bootstrap
            boots{n} = bootseq; 
        end
        
        % Função anônima para construir a árvore usando Neighbor-Joining e calcular a distância
        fun = @(x) seqlinkage(x,'average',{primates.Header});

        % Inicializa células para armazenar as árvores bootstrap
        boot_trees = cell(num_boots,1); 
       
             
        % Calcula a distância entre as replicatas bootstrap e constrói as árvores correspondentes
        for n = 1:num_boots            
            dist_tmp = seqpdist(boots{n});
            boot_trees{n} = fun(dist_tmp);
        end
    
        
        % Obtém a topologia e as espécies para cada ramo da árvore original
        for i = num_seqs-1:-1:1  % for every branch, reverse order to preallocate
            branch_pointer = i + num_seqs;
            sub_tree = subtree(NJtree,branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_species{i} = sort(get(sub_tree,'LeafNames'));
        end
        
        % Obtém a topologia e as espécies para cada ramo das árvores bootstrap
        for j = num_boots:-1:1
            for i = num_seqs-1:-1:1  % for every branch
                branch_ptr = i + num_seqs;
                % deve ser esta linha
                sub_tree = subtree(boot_trees{j},branch_ptr);
                clusters_pointers{i,j} = getcanonical(sub_tree);
                clusters_species{i,j} = sort(get(sub_tree,'LeafNames'));
            end
        end
       
        % Contagem de ramos com topologia similar entre a árvore original e as árvores bootstrap
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

        % Cálculo da probabilidade de confiança (Pc)
        Pc = count ./ num_boots;   % confidence probability (Pc)
    
        % Exibição dos valores de confiança na árvore original
        [ptrs,dist,names] = get(NJtree,'POINTERS','DISTANCES','NODENAMES');
        
        for i = 1:num_seqs -1  % for every branch
            branch_ptr = i + num_seqs;
            names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
        end
    
        % Constrói a árvore final   
        tr = phytree(ptrs,dist,names); 
    

        % Mostrar os pontos de ramificação que têm um nivel de confiança maior que C
        % Um nível de confiança mais alto (por exemplo, 99%) possa fornecer uma estimativa mais precisa, 
        % Um nível de confiança mais baixo (por exemplo, 90%) pode resultar 
        % em intervalos de confiança mais amplos e menos precisos.        
        c = input('nivel de confiança maior que: ');
        
        % Verificando se c está dentro do intervalo
        if c >= 0 && c <= 1
             high_conf_branch_ptr = find(Pc > c) + num_seqs;
             view(tr, high_conf_branch_ptr);
        else
            disp('Insira novamente o nível de confiança, valor inválido.');
        end

        disp('Processo concluído com sucesso.');
    end
end

disp('')