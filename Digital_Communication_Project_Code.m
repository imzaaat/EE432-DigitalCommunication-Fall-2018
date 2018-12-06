clc;
clear all;
set(0,'RecursionLimit',1000)
format long
format compact

% Generate 5 random numbers within range.
src_min = 33;
src_max = 126;
n = 5;
src = fix(src_min + rand(1,n) * (src_max - src_min));

% Generating random probabilities.
p1 = 0.8 + rand(1) * (0.1);
p2 = (0.08 + rand(1) * (0.01));
p3 = (0.008 + rand(1) * (0.001));
p4 = (0.0008 + rand(1) * (0.0001));
ps = [p1,p2,p3,p4];
p5 = 1 - sum(ps);
prob = [p1 p2 p3 p4 p5];

% Sort the probabilities and reorder the source accordingly.
[prob, prob_sorted] = sort(prob, 'descend');
src = src(:,prob_sorted);

% Calculating entropy of the source.
entropy = 0;
for i = 1:5
    entropy = entropy + prob(i) * log2(1/prob(i));
end

% Display information
fprintf('\t\tSource Information\n')
fprintf('RNG\t\t\tASCII\t\tPROBABILITY\t\n')
fprintf('-----------------------------------\n')

for i = 1:5
   fprintf('%s\t\t\t%s\t\t\t%s\n', num2str(src(i)), char(src(i)), num2str(prob(i))); 
end

fprintf('Source Entropy: %f\n', entropy)


% Computing the n-gram for n = 1, 2, 3, and 4.
% Flip function is used to enhance readability of the combination, reverses
% order of elements.

%%%%%%%%%%
fprintf('\n\n\t\tN-GRAM FOR n = 1\n\n\n')
%%%%%%%%%%

n1 = src';
n1p = prob';

[n1_codes, n1_avgln] = huffmandict(n1, n1p);

fprintf('PROB\t\tSYMBOL\t\tCODEWORD\n')
for i = 1:5
   fprintf('%f\t\t%s\t\t\t%s\n', n1p(i), n1(i, :), num2str(cell2mat(n1_codes(i,2))))  
end

fprintf('Average Length: %f (bits/symbol)\n', n1_avgln)
fprintf('Efficiency: %%%f\n', entropy/n1_avgln*100)

%%%%%%%%%%
fprintf('\n\n\t\tN-GRAM FOR n = 2\n\n\n')
%%%%%%%%%%


n2 = flip(combvec(src, src)', 2);
n2p = [];
counter = 1;
for r1 = 1:5
    for r2 = 1:5
        n2p(counter) = prob(r1) * prob(r2);
        counter = counter + 1;
    end
end

[n2p, n2p_order] = sort(n2p, 'descend');
n2 = n2(n2p_order,:);

[n2_codes, n2_avgln] = huffmandict(1:25, n2p);

fprintf('PROB\t\tSYMBOL\t\t\tCODEWORD\n')
for i = 1:25
   fprintf('%f\t\t%s\t\t\t%s\n', n2p(i), n2(i, :), num2str(cell2mat(n2_codes(i,2)))) 
end

fprintf('Average Length: %f (bits/symbol)\n', n2_avgln/2)
fprintf('Efficiency: %%%f\n', entropy/(n2_avgln/2)*100)

%%%%%%%%%%
fprintf('\n\n\t\tN-GRAM FOR n = 3\n\n\n')
%%%%%%%%%%


n3 = flip(combvec(src, src, src)', 2);
n3p = [];
counter = 1;
for r1 = 1:5
    for r2 = 1:5
        for r3 = 1:5
            n3p(counter) = prob(r1) * prob(r2) * prob(r3);
            counter = counter + 1;
        end
    end
end

[n3p, n3p_order] = sort(n3p, 'descend');
n3 = n3(n3p_order,:);

[n3_codes, n3_avgln] = huffmandict(1:125, n3p);

fprintf('PROB\t\tSYMBOL\t\t\tCODEWORD\n')
for i = 1:125
   fprintf('%f\t\t%s\t\t\t%s\n', n3p(i), n3(i, :), num2str(cell2mat(n3_codes(i,2)))) 
end

fprintf('Average Length: %f (bits/symbol)\n', n3_avgln/3)
fprintf('Efficiency: %%%f\n', entropy/(n3_avgln/3)*100)


%%%%%%%%%%
fprintf('\n\n\t\tN-GRAM FOR n = 4\n\n\n')
%%%%%%%%%%

n4 = flip(combvec(src, src, src, src)', 2);

n4p = [];
counter = 1;
for r1 = 1:5
    for r2 = 1:5
        for r3 = 1:5
            for r4 = 1:5
                n4p(counter) = prob(r1) * prob(r2) * prob(r3) * prob(r4);
                counter = counter + 1;
            end
        end
    end
end

[n4p, n4p_order] = sort(n4p, 'descend');
n4 = n4(n4p_order,:);

[n4_codes, n4_avgln] = huffmandict(1:625, n4p);

fprintf('PROB\t\tSYMBOL\t\t\tCODEWORD\n')
for i = 1:625
   fprintf('%f\t\t%s\t\t\t%s\n', n4p(i), n4(i, :), num2str(cell2mat(n4_codes(i,2))));
end

fprintf('Average Length: %f (bits/symbol)\n', n4_avgln/4)
fprintf('Efficiency: %%%f\n', entropy/(n4_avgln/4)*100)

