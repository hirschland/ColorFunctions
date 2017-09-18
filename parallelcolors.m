function [] = parallelcolors(n,I,L)
% Demonstrates the higher dimensionality of physical spectral colours by
% presenting them in parallel coordinate space.

% Each loop (repeated l times) generates n spectral colours from random decimals.
% In each generation, randomly selected wavelengths are boosted I times by
% predetermined factors. Results are plotted on parallel coordinates in
% order of increasing wavelengths (the dimensions of each colour). Boosting
% more wavelengths will desaturate individual colours.

% Joshua Harvey 2017

if nargin==0
    L=10;n=10;I=10;
end

figure
for l = 1:L    
    %V = linspace(380,700,v);
    V = optgetpref('WLRange');
    
    % random allocation of freqs between 0 and 1
    M = rand([n,length(V)]);
    B = [2 4 6 8 10 15 20]; % boost factors
    M=repmat(M,[1 1 length(B)]);
    for b = 1:length(B)
        for m = 1:n %n
            MB(m,:,b) = M(m,:,b);
            for i = 1:I % number of wavelengths for each colour to boost
                r = ceil(rand(1)*31);
                MB(m,r,b) = M(m,r,b)*B(b);
            end
        end
        %MB(MB<1)=0; % put all non-boosted wavelengths to 0
        D{b} = roo2rgb(MB(:,:,b),[],V);
        D{b} = D{b}(:,:)/max(max(D{b}));
        
        P{b} = parallelcoords(MB(:,:,b));
        
        for p = 1:n
            P{b}(p).Color = [D{b}(p,:) 0.15];
            P{b}(p).LineWidth = 3;
        end
        
        hold on
        
    end
    pause(0.1)
end
end