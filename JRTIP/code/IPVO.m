clear all;
clc

addpath(genpath('Origin Images'));
I = double(imread('Lena.bmp'));
[A B] = size(I);

index_final = 100;

MMI_Man = zeros(2,index_final);
R = zeros(5,5,index_final);

for a = 4:4
    for b =3:3
        [a b]
        LM = zeros(floor(A/a),floor(B/b));
        dmax = zeros(floor(A/a),floor(B/b));
        dmin = zeros(floor(A/a),floor(B/b));
        NL = zeros(floor(A/a),floor(B/b));
        for i = 1:floor(A/a)
            for j = 1:floor(B/b)
                X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
                X = X(:);
                [Y In] = sort(X);
                u = min(In(a*b),In(a*b-1));
                v = max(In(a*b),In(a*b-1));
                dmax(i,j) = X(u)-X(v);
                u = min(In(2),In(1));
                v = max(In(2),In(1));
                dmin(i,j) = X(u)-X(v);
                NL(i,j) = Y(a*b-1)-Y(2);
                if Y(a*b) == 255 || Y(1) == 0
                    LM(i,j) = 1;
                end
            end
        end
%         LMC = LM(:);
%         xC = cell(1,1);
%         xC{1} = LMC;
%         data = Arith07(xC);
        NLmax = max(NL(:))
        ECtemp = 0;%12+2*ceil(log2(A/a*B/b))+8*length(data);
        indextag = 0;
        for index = 5:index_final
            if indextag == 0
                EC = 1000*index;
                for t = 1:NLmax+1
                    PSNR = 0;
                    ECbis = EC+ECtemp;
                    for i = 1:floor(A/a)
                        for j = 1:floor(B/b)
                            if ECbis > 0 && NL(i,j) < t && LM(i,j) == 0
                                if dmax(i,j) == 1 || dmax(i,j) == 0
                                    ECbis = ECbis-1;
                                    PSNR = PSNR+0.5;
                                else
                                    PSNR = PSNR+1;
                                end
                                if dmin(i,j) == 1 || dmin(i,j) == 0
                                    ECbis = ECbis-1;
                                    PSNR = PSNR+0.5;
                                else
                                    PSNR = PSNR+1;
                                end
                            end
                        end
                    end
                    if ECbis <= 0
                        R(a,b,index) = 10*log10(A*B*255^2/PSNR);
                        break
                    end
                end
                if ECbis > 0
                    indextag = 1;
                    index
                end
            end
        end
    end
end

for index = 5:index_final
    temp = max(max(R(:,:,index)));
    if temp > 0
        indexend = index;
        MMI_Man(1,index) = index*1000;
        MMI_Man(2,index) = temp;
    end
end

MMI_Man = MMI_Man(:,5:indexend);
res = MMI_Man;
save('IPVO 4x3_2013_lena', 'MMI_Man')