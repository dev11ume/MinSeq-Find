% SCORE SEQUENCES USING MINSEQS WITH GAP SCORING V1

function int_chip_seq=get_seq_scores_gp_v1(seq_colp5n,l,...
    min_half_nmer,max_half_nmer,max_gap,int_nu2x_cs)

% seq_colp5nr=5-seq_colp5n(:,end:-1:1);
% % REMOVE Ns
% seq_colp5nr(seq_colp5nr==5)=0;

% COUNTING half_nmer GAP half_nmer
[l_seq,len]=size(seq_colp5n);
int_chip_seq=zeros(l_seq,1);

for half_nmeri=min_half_nmer:max_half_nmer
    fivesi = power(5,0:half_nmeri-1);
    fivesix = power(5,0:half_nmeri);
    for half_nmerj=min_half_nmer:max_half_nmer
        fprintf('Nmer-Gap-Nmer=%d\t%d\n',half_nmeri,half_nmerj);
        fivesj = power(5,0:half_nmerj-1);
        fivesjx = power(5,0:half_nmerj);
        for j=0:min(l-half_nmeri-half_nmerj,max_gap)
            % fprintf('Nmer-Gap-Nmer=%d\t%d\t%d\n',half_nmeri,j,half_nmerj);
            tot_l=half_nmeri+half_nmerj+j;

            for i=1:len-tot_l+1
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;

                % NOW LOOKING FOR INX OF LONGER SEQUENCE - IF EXISTS
                if i>1 && tot_l<l && half_nmeri<max_half_nmer
                    inx_LL=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i-1:i+half_nmeri-1).*fivesix(ones(l_seq,1),:),2)+...
                        +(5^(min(max_half_nmer,j*1000+(half_nmeri+1))))*sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                        +1;
                    % DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
                else
                    inx_LL=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                if (i<len-tot_l+1) && (tot_l<l) && half_nmerj<max_half_nmer
                    inx_RR=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                        +(5^(min(max_half_nmer,j*1000+half_nmeri)))*sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj).*fivesjx(ones(l_seq,1),:),2)+...
                        +1;
                    % DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
                else
                    inx_RR=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                if j>0 &&  half_nmeri<max_half_nmer
                    inx_LR=(j-1)*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri).*fivesix(ones(l_seq,1),:),2)+...
                        +(5^(min(max_half_nmer,(j-1)*1000+(half_nmeri+1))))*sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                        +1;
                else
                    inx_LR=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                if j>0 &&  half_nmerj<max_half_nmer
                    inx_RL=(j-1)*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                        +(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*sum(seq_colp5n(:,i+half_nmeri+j-1:i+half_nmeri+j+half_nmerj-1).*fivesjx(ones(l_seq,1),:),2)+...
                        +1;
                    % DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
                else
                    inx_RL=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                inxx=(int_nu2x_cs(inx_LL)>0)|(int_nu2x_cs(inx_RR)>0)|...
                    (int_nu2x_cs(inx_LR)>0)|(int_nu2x_cs(inx_RL)>0);
                inx(inxx)=1;% ASSUMING int_nu2x_cs(1)=0;
                int_chip_seq=max(int_chip_seq,int_nu2x_cs(inx));
            end
        end
    end
end


end
