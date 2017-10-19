GB Comments:
Darlanminussi
Prob1: 100%
Prob2:
P1:100
P2: 100
P3:100
P4:100
P5:0 No Answer provided
Prob3
P1: 100
P2:100
P3:100 
Overall: 88


% Homework 1. Due before class on 9/5/17
% Darlan Conterno Minussi

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
% x = 3; y = '5'; %mixed

%your code goes here

if (ischar(x))
   x = str2num(x);
end
if (ischar(y))
   y = str2num(y);
end
    
%output your answer
disp(x + y);

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length

bases = ['A' 'T' 'G' 'C'];
gen_numbers = randi(numel(bases),[1 N]);
rand_seq = bases(gen_numbers);


%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

N = 500; % define sequence length

bases = ['A' 'T' 'G' 'C'];
gen_numbers = randi(numel(bases),[1 N]);
rand_seq = bases(gen_numbers);

start = strfind(rand_seq, 'ATG');
stop = [strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG')];

start = strfind(rand_seq, 'ATG');
stop = [strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG')];

first_stop_codon = zeros(1,length(start));
for i = 1:length(start)
    length_orf = stop - start(i);
    good_length = 1e8;
    index = 0;
    for j = 1:length(length_orf)
        if (length_orf(j) > 0) & (mod(length_orf(j),3) == 0) & (length_orf(j) < good_length)
            good_length = length_orf(j);
            index = j;
        end
    end
    if index > 0
        first_stop_codon(i) = stop(index);
    else
        first_stop_codon(i) = start(i);
    end
    size_orf = first_stop_codon - start + 3;
    [max_length_orf, idx_max] = max(size_orf);
    
    ORFlength = max_length_orf;
    start_pos = start(idx_max);
    stop_pos = first_stop_codon(idx_max);
end

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.


length_orf_prob = 0;

for iteration = 1:1000
    
    
    N = 500; % define sequence length
    
    bases = ['A' 'T' 'G' 'C'];
    gen_numbers = randi(numel(bases),[1 N]);
    rand_seq = bases(gen_numbers);
    
    start = strfind(rand_seq, 'ATG');
    stop = [strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG')];
    
    
    for i = 1:length(start)
        for j = 1:length(stop)
            length_orf(i) = stop(j) - start(i);
        end
    end
    
    if (~isempty(length_orf))
        length_orf = length_orf(length_orf > 50);
        length_orf = length_orf(mod(length_orf,3) == 0);
    else
    end
 
    
if length(length_orf > 0)
    length_orf_prob = length_orf_prob + 1;
end
    
end


prob = length_orf_prob/1000;
disp(prob);


%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 

length_orf_prob = 0
seq_counter = 1
table = []


for n = 50:10:500
    
    bases = ['A' 'T' 'G' 'C'];
    gen_numbers = randi(numel(bases),[1 n]);
    rand_seq = bases(gen_numbers);
    
    
    start = strfind(rand_seq, 'ATG');
    stop = [strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG')];
    
    
    for iteration = 1:1000
       
        
        bases = ['A' 'T' 'G' 'C'];
        gen_numbers = randi(numel(bases),[1 n]);
        rand_seq = bases(gen_numbers);
        
        start = strfind(rand_seq, 'ATG');
        stop = [strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG')];
        
        
        for i = 1:length(start)
            for j = 1:length(stop)
                length_orf(i) = stop(j) - start(i);
            end
        end
        
        if (~isempty(length_orf))
            length_orf = length_orf(length_orf > 50);
            length_orf = length_orf(mod(length_orf,3) == 0);
        else
        end
        
        
        if length(length_orf > 0)
            length_orf_prob = length_orf_prob + 1;
        end
        
    end
    
    
    prob = length_orf_prob/1000;
    
    table(seq_counter,1) = n
    table(seq_counter,2) = prob
    
    length_orf_prob = 0;
    
    seq_counter = seq_counter + 1;
    
end

scatter(table(:,1), table(:,2));
xlabel('Random sequence N');
ylabel('Frequency orf > 50 bp');

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

% specifying full path, please remove path if testing and make sure file is in the same folder)
rawdata = readtable('/Users/dcminussi/Documents/GitHub/hw1-darlanminussi/qPCRdata.txt');
cp = rawdata(:,{'Cp'});
cp.array = table2array(cp);



% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 

cp = vec2mat(cp.array,12);
cp = cp(1:6,:);

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

norm_genes = cp(:,10:12);
norm_genes_mean = mean(norm_genes, 1);
test = [];
final_results = [];

% calculating mean for the triplicates
for ii = 0:3
    cp_mean(:,ii + 1) = mean(cp(:,ii*3+1:(ii+1)*3),2);
end

for i = (1:size(cp_mean,1))
    for j = (1:size(cp_mean,2)-1)
    test(i,j) = cp_mean(1,j) - cp_mean(i,j) - (cp_mean(1,4) - cp_mean(i,4));
    end
end

final_results = 2.^test;
disp(final_results);

bar(final_results);
ylabel('Fold Change');
xlabel('Conditions');

    
%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


