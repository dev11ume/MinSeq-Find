% SCORE SEQUENCES USING MINSEQS WITH GAP SCORING V2

function int_chip_seq=get_seq_scores_gp_v2(seq_colp5n,l,...
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
                    %                     DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
                else
                    inx_LL=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                if (i<len-tot_l+1) && (tot_l<l) && half_nmerj<max_half_nmer
                    inx_RR=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                        +(5^(min(max_half_nmer,j*1000+half_nmeri)))*sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj).*fivesjx(ones(l_seq,1),:),2)+...
                        +1;
                    %                     DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
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
                    %                     DON'T HAVE TO DO REV COMP BCOZ OF int_nu2x_cs ALREADY HAS THAT
                else
                    inx_RL=ones(l_seq,1);% ASSUMING int_nu2x_cs(1)=0;
                end
                inxx=(int_nu2x_cs(inx_LL)>0)|(int_nu2x_cs(inx_RR)>0)|...
                    (int_nu2x_cs(inx_LR)>0)|(int_nu2x_cs(inx_RL)>0);
                inxn1=find(inxx~=1);
                inx_max_fit=find(int_nu2x_cs(inx(inxn1))>0); % THOSE WHICH GOT MAX FIT WITH NONZERO INT
                inx_max_fit=inxn1(inx_max_fit);
                %                 int_chip_seq(inx_max_fit)=max(int_chip_seq(inx_max_fit),int_nu2x_cs(inx(inx_max_fit)));
                l_seq1=length(inx_max_fit);
                % GET INTENSITY FOR inx_max_fit BY SUBTRACTING OUT
                if l_seq1>0
                    if half_nmeri==max_half_nmer && half_nmerj==max_half_nmer
                        % CAN MAKE FASTER BY ACCESSING int_nu2x_cs ONLY ONCE
                        int_chip_seq(inx_max_fit)=max(int_chip_seq(inx_max_fit),int_nu2x_cs(inx(inx_max_fit)));
                    else
                        if i>1 && tot_l<l && half_nmeri<max_half_nmer
                            inxLLt=j*(5^(2*max_half_nmer))+sum(seq_colp5n(inx_max_fit,i:i+half_nmeri-1).*fivesix(ones(l_seq1,1),2:end),2)+...
                                +(5^(min(max_half_nmer,j*1000+(half_nmeri+1))))*sum(seq_colp5n(inx_max_fit,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq1,1),:),2)+...
                                +1;
                            inxLL1=1+inxLLt;
                            inxLL2=2+inxLLt;
                            inxLL3=3+inxLLt;
                            inxLL4=4+inxLLt;
                            inx_longerLL=(int_nu2x_cs(inxLL1)>0)+(int_nu2x_cs(inxLL2)>0)+(int_nu2x_cs(inxLL3)>0)+(int_nu2x_cs(inxLL4)>0);
                        else
                            inx_longerLL=zeros(l_seq1,1);
                        end
                        if (i<len-tot_l+1) && (tot_l<l) && half_nmerj<max_half_nmer
                            inxRRt=j*(5^(2*max_half_nmer))+sum(seq_colp5n(inx_max_fit,i:i+half_nmeri-1).*fivesi(ones(l_seq1,1),:),2)+...
                                +(5^(min(max_half_nmer,j*1000+half_nmeri)))*sum(seq_colp5n(inx_max_fit,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq1,1),:),2)+...
                                +1;
                            inxRR1=inxRRt+(5^(min(max_half_nmer,j*1000+half_nmeri)))*(1*(5^half_nmerj));
                            inxRR2=inxRRt+(5^(min(max_half_nmer,j*1000+half_nmeri)))*(2*(5^half_nmerj));
                            inxRR3=inxRRt+(5^(min(max_half_nmer,j*1000+half_nmeri)))*(3*(5^half_nmerj));
                            inxRR4=inxRRt+(5^(min(max_half_nmer,j*1000+half_nmeri)))*(4*(5^half_nmerj));
                            inx_longerRR=(int_nu2x_cs(inxRR1)>0)+(int_nu2x_cs(inxRR2)>0)+(int_nu2x_cs(inxRR3)>0)+(int_nu2x_cs(inxRR4)>0);
                        else
                            inx_longerRR=zeros(l_seq1,1);
                        end
                        if j>0 &&  half_nmeri<max_half_nmer
                            inxLRt=(j-1)*(5^(2*max_half_nmer))+sum(seq_colp5n(inx_max_fit,i:i+half_nmeri-1).*fivesi(ones(l_seq1,1),:),2)+...
                                +(5^(min(max_half_nmer,(j-1)*1000+(half_nmeri+1))))*sum(seq_colp5n(inx_max_fit,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq1,1),:),2)+...
                                +1;
                            inxLR1=inxLRt+1*(5^half_nmeri);
                            inxLR2=inxLRt+2*(5^half_nmeri);
                            inxLR3=inxLRt+3*(5^half_nmeri);
                            inxLR4=inxLRt+4*(5^half_nmeri);
                            inx_longerLR=(int_nu2x_cs(inxLR1)>0)+(int_nu2x_cs(inxLR2)>0)+(int_nu2x_cs(inxLR3)>0)+(int_nu2x_cs(inxLR4)>0);
                        else
                            inx_longerLR=zeros(l_seq1,1);
                        end
                        if j>0 &&  half_nmerj<max_half_nmer
                            inxRLt=(j-1)*(5^(2*max_half_nmer))+sum(seq_colp5n(inx_max_fit,i:i+half_nmeri-1).*fivesi(ones(l_seq1,1),:),2)+...
                                +(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*sum(seq_colp5n(inx_max_fit,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesjx(ones(l_seq1,1),2:end),2)+...
                                +1;
                            inxRL1=inxRLt+(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*(1);
                            inxRL2=inxRLt+(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*(2);
                            inxRL3=inxRLt+(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*(3);
                            inxRL4=inxRLt+(5^(min(max_half_nmer,(j-1)*1000+half_nmeri)))*(4);
                            inx_longerRL=(int_nu2x_cs(inxRL1)>0)+(int_nu2x_cs(inxRL2)>0)+(int_nu2x_cs(inxRL3)>0)+(int_nu2x_cs(inxRL4)>0);
                        else
                            inx_longerRL=zeros(l_seq1,1);
                        end

                        inx_nolong=find(inx_longerLL+inx_longerRR+inx_longerLR+inx_longerRL==0);
                        int_chip_seq(inx_max_fit(inx_nolong))=max(int_chip_seq(inx_max_fit(inx_nolong)),int_nu2x_cs(inx(inx_max_fit(inx_nolong))));
                        inx_long2=find(inx_longerLL|inx_longerRR|inx_longerLR|inx_longerRL>0);
                        if ~isempty(inx_long2) % IF THERE EXIST INTENSITY FOR LONGER SEQ NOT MATCHING seq_colp5n
                            intx=int_nu2x_cs(inx(inx_max_fit(inx_long2)));
                            intxFIX=intx;
                            inx_longerLLix=find(inx_longerLL(inx_long2));
                            inx_longerLLix2=inx_long2(inx_longerLLix);
                            if ~isempty(inx_longerLLix2)
                                int_totLL=int_nu2x_cs(inxLL1(inx_longerLLix2))+int_nu2x_cs(inxLL2(inx_longerLLix2))+int_nu2x_cs(inxLL3(inx_longerLLix2))+int_nu2x_cs(inxLL4(inx_longerLLix2));
                                int_LL=(4*intxFIX(inx_longerLLix)-int_totLL)./(4-inx_longerLL(inx_longerLLix2));
                                intx(inx_longerLLix)=min(intx(inx_longerLLix),int_LL);
                            end
                            inx_longerRRix=find(inx_longerRR(inx_long2));
                            inx_longerRRix2=inx_long2(inx_longerRRix);
                            if ~isempty(inx_longerRRix2)
                                int_totRR=int_nu2x_cs(inxRR1(inx_longerRRix2))+int_nu2x_cs(inxRR2(inx_longerRRix2))+int_nu2x_cs(inxRR3(inx_longerRRix2))+int_nu2x_cs(inxRR4(inx_longerRRix2));
                                int_RR=(4*intxFIX(inx_longerRRix)-int_totRR)./(4-inx_longerRR(inx_longerRRix2));
                                intx(inx_longerRRix)=min(intx(inx_longerRRix),int_RR);
                            end
                            inx_longerLRix=find(inx_longerLR(inx_long2));
                            inx_longerLRix2=inx_long2(inx_longerLRix);
                            if ~isempty(inx_longerLRix2)
                                int_totLR=int_nu2x_cs(inxLR1(inx_longerLRix2))+int_nu2x_cs(inxLR2(inx_longerLRix2))+int_nu2x_cs(inxLR3(inx_longerLRix2))+int_nu2x_cs(inxLR4(inx_longerLRix2));
                                int_LR=(4*intxFIX(inx_longerLRix)-int_totLR)./(4-inx_longerLR(inx_longerLRix2));
                                intx(inx_longerLRix)=min(intx(inx_longerLRix),int_LR);
                            end
                            inx_longerRLix=find(inx_longerRL(inx_long2));
                            inx_longerRLix2=inx_long2(inx_longerRLix);
                            if ~isempty(inx_longerRLix2)
                                int_totRL=int_nu2x_cs(inxRL1(inx_longerRLix2))+int_nu2x_cs(inxRL2(inx_longerRLix2))+int_nu2x_cs(inxRL3(inx_longerRLix2))+int_nu2x_cs(inxRL4(inx_longerRLix2));
                                int_RL=(4*intxFIX(inx_longerRLix)-int_totRL)./(4-inx_longerRL(inx_longerRLix2));
                                intx(inx_longerRLix)=min(intx(inx_longerRLix),int_RL);
                            end
                            int_chip_seq(inx_max_fit(inx_long2))=max(int_chip_seq(inx_max_fit(inx_long2)),intx);
                        end
                    end
                end
            end
        end
    end
end


end
