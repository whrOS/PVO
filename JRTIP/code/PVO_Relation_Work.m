%%%%%%%% PVO的改进版本，Xiaolong Li，2018-07-22
% 按照T的增长，没有考虑具体的嵌入容量

tic
clear all;
clc

addpath(genpath('Origin Images'));

I = double(imread('Lena.bmp'));
[A B] = size(I);

index = 0;
R = zeros(7,10^7);

hist2D = zeros(512,512);

hist1D_dmax = zeros(1,512);
hist1D_dmin = zeros(1,512);

for a = 3:3
    for b = 3:3
        [a b];
        NL = zeros(floor((A-2)/a),floor((B-2)/b));
        %         EC = zeros(floor((A-2)/a),floor((B-2)/b));
        EC=0;
        %         ED = zeros(floor((A-2)/a),floor((B-2)/b));
        ED=0;
        for i = 1:floor((A-2)/a)
            for j = 1:floor((B-2)/b)
                for ii = 1:a+1
                    for jj = 1:b+2
                        if ii == a+1 || jj == b+1 || jj == b+2
                            NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii+1,b*(j-1)+jj));
                        end
                    end
                end
                for ii = 1:a+2
                    for jj = 1:b+1
                        if ii == a+1 || ii == a+2 || jj == b+1
                            NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii,b*(j-1)+jj+1));
                        end
                    end
                end
                X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
                X = X(:);
                [Y In] = sort(X);
                % max
                dmax = Y(a*b) - Y(a*b-1);
                % min
                dmin = Y(1) - Y(2);
                hist1D_dmax(dmax+1) = hist1D_dmax(dmax+1) + 1;
                hist1D_dmin(dmin+256) = hist1D_dmin(dmin+256) + 1;
            end
            
        end
    end
end
hist1D_dmax=hist1D_dmax(1:max(find(hist1D_dmax>=50)));
hist1D_dmin=hist1D_dmin(min(find(hist1D_dmin>=20)):max(find(hist1D_dmin>=20)));
bar(0:19,hist1D_dmax);
axis([-1, 20, 0, 8100]);
