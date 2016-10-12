function H_dct = mdct_filterbank(M)

    h = [];

    vn=0:(M*2-1);
    w = sin((pi/(4*M))*(2*vn+1)); % window function for perfect reconstruction.

    H_dct=[];

    for k=0:M-1
        
        for n=0:M*2-1
            
            H_dct(k+1,n+1)=w(n+1) * sqrt(2/M) * cos( (2*n+M+1)*(2*k+1)*pi / (4*M) );
            
        end
    end
