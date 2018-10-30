%%%%%%%% PVO�ĸĽ�汾��Xiaolong Li��2018-07-22
% ����T��������û�п��Ǿ����Ƕ������?

tic
clear all;
clc

addpath(genpath('Origin Images')); 

Iname = 'Airplane';
I = double(imread([Iname '.bmp']));
[A B] = size(I);

index = 0;
R = zeros(7,10^7);

hist2D = zeros(512,512);

T = 20;%% 5000bit WPVO 62, IPVO 29
hist2D_flat = zeros(512,512);
hist1D_dmax = zeros(1,512);
hist1D_dmin = zeros(1,512);
hist1D_dmax_all = zeros(1,512);
hist1D_dmin_all = zeros(1,512);

for a = 2:2
    for b = 2:2
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
                if In(a*b) < In(a*b-1)
                    dddd = Y(a*b) - Y(a*b-1) -1;
                    dmax = Y(a*b-1)  - Y(a*b);
                else
                    dddd = Y(a*b) - Y(a*b-1);
                    dmax = Y(a*b) - Y(a*b-1);
                end
                % min
                if In(2) < In(1)
                    ddd = Y(2) - Y(1) -1;
                    dmin = Y(1) - Y(2);
                else
                    ddd = Y(2) - Y(1);
                    dmin = Y(2) - Y(1);
                end
%                 [dmin,dmax]
                hist2D(dmax+256,dmin+256) = hist2D(dmax+256,dmin+256)+1;
                hist1D_dmax_all(dmax+256) = hist1D_dmax_all(dmax+256) + 1;
                hist1D_dmin_all(dmin+256) = hist1D_dmin_all(dmin+256) + 1;
                if NL(i,j) < T
                    if dddd == 0 && ddd == 0
                        EC = EC+log2(3);
                        ED = ED+2/3;
                    else if dddd == 0 || ddd == 0
                           EC = EC+1;
                            ED = ED+1.5;
                        else if dddd == 1 && ddd == 1
                                EC = EC+1;
                                ED = ED+1;
                            else
                                ED = ED+2;
                            end
                        end
                    end
                    hist2D_flat(dmax+256,dmin+256) = hist2D_flat(dmax+256,dmin+256)+1;
                    hist1D_dmax(dmax+256) = hist1D_dmax(dmax+256) + 1;
                    hist1D_dmin(dmin+256) = hist1D_dmin(dmin+256) + 1;
                end

            end
        end
    end
end
idx=min(find(hist1D_dmax_all>=50));
hist1D_dmax_all=hist1D_dmax_all(min(find(hist1D_dmax_all>=50)):max(find(hist1D_dmax_all>=50)));
% hist1D_dmin=hist1D_dmin(min(find(hist1D_dmin>=20)):max(find(hist1D_dmin>=20)));
bar(-(256-idx):1:-(256-idx)+numel(hist1D_dmax_all)-1,hist1D_dmax_all);
axis([-(256-idx)-0.5,-(256-idx)+numel(hist1D_dmax_all)-1+0.5, 0, 6000]);

binS = 5;
x = 256-binS:1:256+binS;
y = 256-binS:1:256+binS;
xbins = -binS:binS;
ybins = -binS:binS;
% hist2D
figure;
% subplot(1,2,1);
bar3(hist2D(x,y),1);
set(gca,'XTickLabel',xbins);
set(gca,'YTickLabel',ybins);
% title(Iname);
% subplot(1,2,2);
% bar3(hist2D_flat(x,y),1);
% set(gca,'XTickLabel',xbins);
% set(gca,'YTickLabel',ybins);

hist2D_ratio = 0;
binS = 5;
for i = 256-binS:1:256+binS
    for j = 256-binS:1:256+binS
        hist2D_ratio = hist2D_ratio + hist2D(i,j);
    end
end
hist2D_ratio = hist2D_ratio / sum(hist2D(:))

% hist2D_flat_ratio = 0;
% binS = 3;
% for i = 256-binS:1:256+binS
%     for j = 256-binS:1:256+binS
%         hist2D_flat_ratio = hist2D_flat_ratio + hist2D_flat(i,j);
%     end
% end
% hist2D_flat_ratio = hist2D_flat_ratio / sum(hist2D_flat(:));

% figure;
% % binS = 10;
% x = 256-binS:1:256+binS;
% y = 256-binS:1:256+binS;
% xbins = -binS:binS;
% ybins = -binS:binS;
% imagesc(hist2D(x,y))
% set(gca,'XTickLabel',xbins);
% set(gca,'YTickLabel',ybins);

% figure;
% % binS = 10;
% x = 256-binS:1:256+binS;
% y = 256-binS:1:256+binS;
% xbins = -binS:binS;
% ybins = -binS:binS;
% 
% subplot(1,2,1);
% bar(xbins,hist1D_dmax(x),1);
% 
% subplot(1,2,2);
% bar(xbins,hist1D_dmin(x),1);
% 

toc

% H = hist2D_flat(x,y);
% EC_00=(H(binS+1,binS+1)+H(binS+1,binS)+H(binS+2,binS)+H(binS+2,binS+1))*log2(3);
% EC_C = sum(sum(H(1:binS,binS:binS+1)))+sum(sum(H(binS+3:end,binS:binS+1))); % 00 0-1��
% EC_R = sum(sum(H(binS:binS+1,1:binS-1)))+sum(sum(H(binS:binS+1,binS+2:end))); % 00��
% EC_Diag = H(binS,binS+2)+H(binS+3,binS+2)+H(binS+3,binS-1)+H(binS,binS-1);
% WPVO_EC = ED_00+ED_C+ED_R+EC_Diag;

% WPVO_ED = (H(binS+1,binS+1)+H(binS+1,binS)+H(binS+2,binS)+H(binS+2,binS+1)) * 2/3 + ...
%     ED_C*1.5+ED_R*1.5+ED_Diag*2;
% HS_Pixels = sum(sum(H)) - ...
%     (H(binS+1,binS+1)+H(binS+1,binS)+H(binS+2,binS)+H(binS+2,binS+1)) - ...
%     sum(sum(H(1:binS,binS:binS+1)))-sum(sum(H(binS+3:end,binS:binS+1))) - ...
%     sum(sum(H(binS:binS+1,1:binS-1)))-sum(sum(H(binS:binS+1,binS+2:end))) - ...
%     (H(binS,binS+2)+H(binS+3,binS+2)-H(binS+3,binS-1)+H(binS,binS-1));
% WPVO_ED = WPVO_ED + HS_Pixels*2;
% 
% IPVO_EC = hist1D_dmax(256)+hist1D_dmax(257)+hist1D_dmin(256)+hist1D_dmin(257);
% IPVO_ED = (hist1D_dmax(256)+hist1D_dmax(257)+hist1D_dmin(256)+hist1D_dmin(257))*0.5 + ...
%     sum(hist1D_dmax(1:255))+sum(hist1D_dmax(258:end))+sum(hist1D_dmin(1:255))+sum(hist1D_dmin(258:end));
% 
% [WPVO_EC WPVO_ED 10*log10(512*512*255^2/WPVO_ED);
%     IPVO_EC IPVO_ED 10*log10(512*512*255^2/IPVO_ED)]