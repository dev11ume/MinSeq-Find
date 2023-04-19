% SCORE SEQUENCES USING MINSEQS WITH NO GAP SCORING V1
function int_chip_seq=get_seq_scores_ng_v1(seq_colp5n,l,...
    min_nmer,max_nmer,int_nu2x_cs)

% seq_colp5nr=5-seq_colp5n(:,end:-1:1);
% % REMOVE Ns
% seq_colp5nr(seq_colp5nr==5)=0;

% COUNTING half_nmer GAP half_nmer
[l_seq,len]=size(seq_colp5n);
int_chip_seq=zeros(l_seq,1);

for nmer=min_nmer:max_nmer
    fprintf('nmer:%d\n',nmer);
    fivesi = power(5,0:nmer-1);
    fivesix = power(5,0:nmer);
    tot_l=nmer;
    for i=1:len-tot_l+1
        inx=sum(seq_colp5n(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;

        % NOW LOOKING FOR INX OF LONGER SEQUENCE - IF EXISTS
        if i>1 && tot_l<l && nmer<max_nmer
            inx_L=sum(seq_colp5n(:,i-1:i+nmer-1).*fivesix(ones(l_seq,1),:),2)+1;
            % DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
        else
            inx_L=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
        end

        if (i<len-tot_l+1) && (tot_l<l) && nmer<max_nmer
            inx_R=sum(seq_colp5n(:,i:i+nmer).*fivesix(ones(l_seq,1),:),2)+1;
            % DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
        else
            inx_R=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
        end

        inxx=(int_nu2x_cs(inx_L)>0)|(int_nu2x_cs(inx_R)>0);
        inx(inxx)=1;% ASSUMING int_nu2x_cs(1)=0;
        int_chip_seq=max(int_chip_seq,int_nu2x_cs(inx));

    end

end


end
