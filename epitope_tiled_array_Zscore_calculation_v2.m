clear
clc
read_singleSample=false(1); % turn it true when Zscores will be used
sampleNo=1; % name of the sample folder containing Z score fits
aminos='ADEFGHKLNPQRSVWY';
alpha=0.05;% this is needed to calculate bonferroni correction
Logfit=true(1);
R_prot=10; %this is the length of sequences that an epitope is cut up into
save_values=true(1);
title_name = 'Prediction of epitopes (IEDB) by Ind ML models';
xlabel_name = 'Chagas proteome sequence';
Proteme_length = 9497648; % length of a proteome
Proteome_tile = true(1);
Find_epitope_Proteome =false(1);
%%

[seqfile,seqpath]=uigetfile('*','Multiselect','off'); % load the epitope list downloaded from IEDB
fid1 = fopen(fullfile([seqpath,seqfile])); 
seq_cell=textscan(fid1,'%s %d %d %s %s','Delimiter',','); % read the list
epitope_raw= seq_cell{1,1}; % 
epitope_raw_char=char(epitope_raw);
start_position = seq_cell{1,2};
end_position = seq_cell{1,3};
antigen_name= seq_cell{1,4};
antigen_name = char(antigen_name);
pathogen_name = seq_cell{1,5};
pathogen_name = char(pathogen_name);

[N,~] = size(epitope_raw);
[tiled_array,seqindex,start_end_Pos]=epitope2peparray(epitope_raw,R_prot,start_position,end_position);
unique_tiledEpitope=unique(cellstr(tiled_array),'stable');
num_unique_tiledEpitope =length(unique_tiledEpitope); % find the number of unique tiled epitopes

% store the tiled array peptide sequence
real_sequence_save = tiled_array;

% now predict the binding of the tiled epitope array
%for the parfor loop, you need dummy versions of these even if we dont use
%them

if ~read_singleSample
    for type_dataset=1:2
        disp('Select the parameter file associated with the set of runs you want to use')
        [Wfile,pathname]=uigetfile('*.txt','MultiSelect','off'); %gives the user a gui to open the file
        param=readtable([pathname,Wfile]);
        param=table2cell(param);
        aminoEigen=str2double(param(find(strcmp(param,'aminoEigen')),2));
        dataShift=str2double(param(find(strcmp(param,'dataShift')),2));
        HiddenLayer=str2double(param(find(strcmp(param,'hiddenLayers')),2));
        aminoEigen=str2double(param(strcmp(param,'aminoEigen'),2));
        % extract dataShift or noise from the parameter
       dataShift= param(find(strcmp(param,'dataShift')),2);
       if strcmp(dataShift{1},cellstr('True'))
          dataShift= true(1);
       elseif strcmp(dataShift{1},cellstr('False'))
           dataShift = false(1); 
       else
           dataShift=str2double(dataShift);
       end
       if type_dataset==1
            HiddenLayer_case=HiddenLayer;
            aminoEigen_case = aminoEigen;
            dataShift_case = dataShift;
            fprintf('AminoEigen used in the NN program is %.f\n',aminoEigen_case)
            fprintf('No of Hidden layers is %d\n',HiddenLayer_case)
       else
            HiddenLayer_control=HiddenLayer;
            aminoEigen_control = aminoEigen;
            dataShift_control = dataShift;
            fprintf('AminoEigen used in the NN program is %.f\n',aminoEigen_control)
            fprintf('No of Hidden layers is %d\n',HiddenLayer_control)
       end 
    
    
    
       % print the filename that was used in the NN program
       filename_index=find(pathname=='\');
       filename=pathname(filename_index(8)+1:filename_index(9)-1);
       fprintf('reading folder %s\n',filename);

       foldername='Sample';
       numsamples=0;

       disp('Reading files');
       %determine how many folders there are
       while exist([pathname,foldername,num2str(numsamples+1)],'dir')==7 %as long as the next folder is there
             numsamples=numsamples+1;
       end


       %read the fit matrices
       disp('reading weights of all samples')
       for isample=1:numsamples
           path=[pathname,foldername,num2str(isample),'\'];
           W(isample)=read_weights_w_bias(path);
       end
       if type_dataset==1
          W_case = W;
          numsample_case = numsamples;
       else
          W_control = W;
          numsample_control = numsamples;
       end
       %Generate an array of sequences to test
       % read data
       fprintf('reading fit data\n');
       [M,K]=size(W(1).W1);
       [RR,H]=size(W(1).W2); % it calculates the length of the peptides from the proteome
       fitR=RR/K;
       % *******Robayet's code ends here.***************

       %***************************************************
       sequence = real_sequence_save;
       [num_seq,~]=size(sequence);
       seq=char(zeros(num_seq,fitR));
       seq(:,1:R_prot)=sequence;
       sequence=seq;
       [N,R]=size(sequence);
       Ntot=N;
       %substitute equivalent amino acids %read in fitting data (I=V, T=S,
       %M=L, C is ignored and not given binding significance).  You have to put
       %something in for these amino acids in a real sequence.  These
       %substitutions are well supported by similarity matrices except for C to A,
       %which shows only weak similarity
       %sequence=real_sequence_save; %store the real sequence
       sequence(sequence=='I')='V';
       sequence(sequence=='T')='S';
       sequence(sequence=='M')='L';
       sequence(sequence=='C')='A';

       seqtot=sequence;
       %project the fits onto a protein
       fprintf('number of peptides used in fit = %d\n',N);
       disp('projecting onto a protein');
  
       F_calc=zeros(N,numsamples);
       parfor isample=1:numsamples
            F_calc(:,isample)=MLevaluate_N_layer_w_bias(sequence,aminos,W(isample));
       end
       delete(gcp('nocreate')) % shut down parallel pool

       if Logfit
           F_calc =10.^(F_calc);
       end
 
        
       if type_dataset==1
           F_calc_case = F_calc;
       else
           F_calc_control = F_calc;
       end
       F_calc = [];
    end
else
    disp('Select the parameter file associated with the set of runs you want to use')
    [Wfile,pathname]=uigetfile('*.txt','MultiSelect','off'); %gives the user a gui to open the file
    param=readtable([pathname,Wfile]);
    param=table2cell(param);
    aminoEigen=str2double(param(find(strcmp(param,'aminoEigen')),2));
    dataShift=str2double(param(find(strcmp(param,'dataShift')),2));
    HiddenLayer=str2double(param(find(strcmp(param,'hiddenLayers')),2));
    fprintf('AminoEigen used in the NN program is %.f\n',aminoEigen)
    fprintf('No of Hidden layers is %d\n',HiddenLayer)
    aminoEigen=str2double(param(find(strcmp(param,'aminoEigen')),2));
    % extract dataShift or noise from the parameter
    dataShift= param(find(strcmp(param,'dataShift')),2);
    if strcmp(dataShift{1},cellstr('True'))
        dataShift= true(1);
    elseif strcmp(dataShift{1},cellstr('False'))
        dataShift = false(1); 
    else
        dataShift=str2double(dataShift);
    end
    fprintf('AminoEigen used in the NN program is %.f\n',aminoEigen)
    fprintf('No of Hidden layers is %d\n',HiddenLayer)

    % print the filename that was used in the NN program
    filename_index=find(pathname=='\');
    filename=pathname(filename_index(8)+1:filename_index(9)-1);
    fprintf('reading folder %s\n',filename);

    foldername='Sample';
    numsamples=1;
    disp('Reading files');
    %read the fit matrices
    fprintf('reding weights of sample no %d\n',sampleNo)
    path=[pathname,foldername,num2str(sampleNo),'\'];
    W=read_weights_w_bias(path);
   %Generate an array of sequences to test
   % read data
   fprintf('reading fit data\n');
   [M,K]=size(W(1).W1);
   [RR,H]=size(W(1).W2); % it calculates the length of the peptides from the proteome
   fitR=RR/K;
   % *******Robayet's code ends here.***************
   %read in the fastA protein or proteome.  You can download these files very
   %simply from uniprot.org.  If you go to a proteome, you can just look for a
   %specific protein, download a whole proteome or select parts.  Note that if
   %you download the whole human proteome and try to calculate C we will run
   %out of memory.  So we will have to break the large proteomes into pieces.
       %***************************************************
       sequence = real_sequence_save;
       [num_seq,~]=size(sequence);
       seq=char(zeros(num_seq,fitR));
       seq(:,1:R_prot)=sequence;
       sequence=seq;
       [N,R]=size(sequence);
       Ntot=N;
       %substitute equivalent amino acids %read in fitting data (I=V, T=S, 
       %M=L, C is ignored and not given binding significance).  You have to put
       %something in for these amino acids in a real sequence.  These
       %substitutions are well supported by similarity matrices except for C to A,
       %which shows only weak similarity
%        real_sequence=sequence; %store the real sequence
       sequence(sequence=='I')='V';
       sequence(sequence=='T')='S';
       sequence(sequence=='M')='L';
       sequence(sequence=='C')='A';

       seqtot=sequence;
       %project the fits onto a protein
       fprintf('number of peptides used in fit = %d\n',N);
       disp('projecting onto a protein');
       Z_calc=MLevaluate_N_layer_w_bias(sequence,aminos,W);
end
    

% calculate Z scores, p_values and Average ratio from two chorts
if ~read_singleSample
    Z_calc=(mean(F_calc_case,2)-mean(F_calc_control,2))./(sqrt(std(F_calc_case,0,2).^2+std(F_calc_control,0,2).^2));

    % calculate p value
    [~,pvalue]=ttest2(F_calc_case',F_calc_control');
    bf_cutoff = alpha/numel(pvalue);
    % calculate average ratio
    ave_ratio=mean(F_calc_case,2)./mean(F_calc_control,2);
end

%print out the best binders
[Z_tiled_Epitope_sort,Zindex2]=sort(Z_calc,'descend');
top_seq=real_sequence_save(Zindex2,:);
top_start_end_pos=start_end_Pos(Zindex2,:);

fprintf('\n*******The highest differential tiled peptides from the IEDB epitope *******\n');
for i=1:10
     fprintf('%s     %.f  %.f  %4.3f  %s\n',top_seq(i,:),top_start_end_pos(i,1),...
     top_start_end_pos(i,2), Z_tiled_Epitope_sort(i), antigen_name(seqindex(Zindex2(i)),:));
end


%plot the Z-scors against the starting position of the tiled array
% peptides
figure(1)
scatter(start_end_Pos(:,1),Z_calc,5)
xlabel(xlabel_name)
ylabel('Z-scores')
xlim([1 Proteme_length]) % change the maximum value accoding to the length of the proteome
title(title_name)

%This is for HCV only
Z_min_max = zeros(length(epitope_raw),2);
Pval_epit= zeros(length(epitope_raw),2);

for jj= 1: length(epitope_raw)
    temp_ind = seqindex==jj; 
    temp_val1_Z = Z_calc(temp_ind);
    Z_min_max(jj,1) = min(temp_val1_Z); 
    Z_min_max(jj,2) = max(temp_val1_Z);
    if ~read_singleSample
        temp_val1_pval = pvalue(temp_ind);
        Pval_epit(jj,1) = min(temp_val1_pval);
        Pval_epit(jj,2) = max(temp_val1_pval);
    end
end

if save_values
    N= length(epitope_raw);
    [savefile,savedir]=uiputfile('*.csv','Input the file name for the sequences, pvalues, average ratios, and Z score');
    disp('Saving average ratio sequences and Z score data in one file')
    temp_file=cell(N+1,1);
for i=1:N
    if ~read_singleSample
        temp_file{1}=sprintf('Epitope,starting position,end position,pvalue_min,pvalue_max,Bf_cutoff,Zscore_min,Zscore_max,Antigen, Pathogen');
        temp_file{i+1}=sprintf('%s,%d,%d,%d,%d,%d,%d,%d,%s,%s',epitope_raw_char(i,:),start_position(i),end_position(i),Pval_epit(i,1),Pval_epit(i,2),bf_cutoff,Z_min_max(i,1),Z_min_max(i,2),antigen_name(i,:),pathogen_name(i,:));
    else
        temp_file{1}=sprintf('Epitope,starting position,end position,Zscore_min,Zscore_max,Antigen, Pathogen');
        temp_file{i+1}=sprintf('%s,%d,%d,%d,%d,%s,%s',epitope_raw_char(i,:),start_position(i),end_position(i),Z_min_max(i,1),Z_min_max(i,2),antigen_name(i,:),pathogen_name(i,:));
    
    end
end

fid5=fopen([savedir,savefile],'w');
for ii=1:N+1
fprintf(fid5, '%s\n', temp_file{ii});
end
fclose(fid5);
clear temp_file
end

    


% plot(start_position(tiled_seqindex_final),Z_final)
%% Save results

% if save_values
%     N= length(Z_calc);
%     [savefile,savedir]=uiputfile('*.csv','Input the file name for the sequences, pvalues, average ratios, and Z score');
%     disp('Saving average ratio sequences and Z score data in one file')
%     temp_file=cell(N+1,1);
% for i=1:N
%     if ~read_singleSample
%         temp_file{1}=sprintf('Sequence,starting position,end position,pvalue,Bf_cutoff,Average_ratio,Zscore,Epitope no,Epitope, Antigen, Pathogen');
%         temp_file{i+1}=sprintf('%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s',real_sequence_save(i,:),start_end_Pos(i,1),start_end_Pos(i,2),pvalue(i),bf_cutoff,ave_ratio(i),Z_calc(i),seqindex(i),epitope_raw_char(seqindex(i),:),antigen_name(seqindex(i),:),pathogen_name(seqindex(i),:));
%     else
%         temp_file{1}=sprintf('Sequence,starting position,end position,Zscore,Epitope no,Epitope, Antigen, Pathogen');
%         temp_file{i+1}=sprintf('%s,%d,%d,%d,%d,%s,%s,%s',real_sequence_save(i,:),start_end_Pos(i,1),start_end_Pos(i,2),Z_calc(i),seqindex(i),epitope_raw_char(seqindex(i),:),antigen_name(seqindex(i),:),pathogen_name(seqindex(i),:));
%         %temp_file{i+1}=sprintf('%s,%d,%d,%d',real_sequence(i,:),start_position(tiled_seqindex_final(i)),end_position(tiled_seqindex_final(i)),Z_final(i));
%     
%     end
% end
% 
% fid5=fopen([savedir,savefile],'w');
% for ii=1:N+1
% fprintf(fid5, '%s\n', temp_file{ii});
% end
% fclose(fid5);
% clear temp_file
% end

%% match the IEDB epitopes with the predicted epitopes by the ML model/s on a proteome of interest

if Find_epitope_Proteome
    
 disp('loading Prediction file') % load the predicted binding data or Z-scores with the tiled peoteome array
[Projection,path1]=uigetfile('*','Multiselect','off');
 fprintf('reading Prediction file %s\n',Projection);
 Projection_data=readtable([path1,Projection]);
 Proteome_seq=Projection_data.Sequence;
 Z_Proteome=Projection_data.Zscore;
 top_predicted_index = Projection_data.PeptidePosition;
 Protein=Projection_data.Protein;


% [Z_Proteome,Zindex_Proteome]=sort(Z_Proteome,'descend');
% top_predicted_index =top_predicted_index(Zindex_Proteome);
% Protein = Protein(Zindex_Proteome);

tiled_array_save= tiled_array;
tiled_array= char(unique_tiledEpitope);
% now locate the known epitope region among the predicted epitopes
match_Protpep_index = false(length(Z_Proteome), length(tiled_array));
pathogen_interest = false(length(Z_Proteome), length(tiled_array));
for ii = 1: length(tiled_array)
    index_match= contains(Proteome_seq, cellstr(tiled_array(ii,:)));
    match_Protpep_index(:,ii)= index_match;
    pathogen_interest(:,ii)= index_match;
end
% find the number, sequences, and ranking of known epitopes and  on the proteome identified the ML model(s)
temp_match = sum(match_Protpep_index,2);
ProtEpit_rank_ML = find(temp_match); % find the ranking of the peptides containing known epitopes on the proteome
total_match = sum(temp_match); % total known epitopes identified by the model
seq_epit_ML = char(Proteome_seq(ProtEpit_rank_ML));% find the sequences
seq_epit_ML_Z = Z_Proteome(ProtEpit_rank_ML); % find the Z-scores
seq_epit_ML_position =top_predicted_index(ProtEpit_rank_ML); % find the positions on the proteome
seq_epit_ML_Protein = char(Protein(ProtEpit_rank_ML));
%save these information 
if save_values
    Nepit= length(ProtEpit_rank_ML );
    [savefile,savedir]=uiputfile('*.csv','Input the file name for the sequences, pvalues, average ratios, and Z score');
    disp('Saving information of identified epitopes on the proteomes in one file')
    temp_file=cell(N+1,1);
for i=1: Nepit
        temp_file{1}=sprintf('Sequence,Proteome_starting position,Z_rank,Zscore,Protein');
        temp_file{i+1}=sprintf('%s,%d,%d,%d,%s',seq_epit_ML(i,:),seq_epit_ML_position(i),ProtEpit_rank_ML(i),seq_epit_ML_Z(i),seq_epit_ML_Protein(i,:));       
end

fid6=fopen([savedir,savefile],'w');
for ii=1:N+1
fprintf(fid6, '%s\n', temp_file{ii});
end
fclose(fid6);
clear temp_file
end



% calculate enrichment scores using the total number of available epitopes
% and the predicted known epitopes identified in the proteome by the ML
% model(s)
% num_unique_tiledEpitope =length(unique(cellstr(tiled_array))); % find the number of unique tiled epitopes
% fraction_IEDB_Proteome = num_unique_tiledEpitope/length(Proteome_seq); % enrichment using known epitopes
% prob_epitope_MLpred = total_match/length(Proteome_seq); % enrichment using known predicted epitopes
% 
% % visulalize the predicted z-scores on the epitopes on the proteome
% figure(2)
% [position_sort, tempx]= sort(top_predicted_index);
% Z_prot_temp = Z_Proteome(tempx);
% scatter(ProtEpit_rank_ML,Z_prot_temp(ProtEpit_rank_ML),5);
% xlabel(xlabel_name)
% ylabel('Z-scores')
% xlim([1 Proteme_length]) % change the maximum value accoding to the length of the proteome
% title(title_name)



% this will calculate how many of the unique tiled epitopes are part of the
% proteome used for projection 
temp = false(length(tiled_array),1);
for i =1:length(tiled_array)
temp_index=contains(Proteome_seq,cellstr(tiled_array(i,:)));
if sum(temp_index)>0
temp(i)=true(1);
end
end

% find the number of unique epitopes in the proteome
Uniq_epitope= unique(cellstr(tiled_array(temp,:)),'stable');
num_uniq_epitope = length(Uniq_epitope);
end



