% SCORE SEQUENCES USING MINSEQS WITH NO GAP SCORING V2

function int_chip_seq=get_seq_scores_ng_v2(seq_colp5n,l,...
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

        inxn1=find(inxx~=1);
        inx_max_fit=find(int_nu2x_cs(inx(inxn1))>0); % THOSE WHICH GOT MAX FIT WITH NONZERO INT
        inx_max_fit=inxn1(inx_max_fit);
        l_seq1=length(inx_max_fit);
        % GET INTENSITY FOR inx_max_fit BY SUBTRACTING OUT
        if l_seq1>0
            if nmer==max_nmer
                % CAN MAKE FASTER BY ACCESSING int_nu2x_cs ONLY ONCE
                int_chip_seq(inx_max_fit)=max(int_chip_seq(inx_max_fit),int_nu2x_cs(inx(inx_max_fit)));
            else
                if i>1 && tot_l<l && half_nmeri<max_half_nmer
                    inxLLt=sum(seq_colp5n(inx_max_fit,i-1:i+nmer-1).*fivesix(ones(l_seq1,1),2:end),2)+1;
                    inxLL1=1+inxLLt;
                    inxLL2=2+inxLLt;
                    inxLL3=3+inxLLt;
                    inxLL4=4+inxLLt;
                    inx_longerLL=(int_nu2x_cs(inxLL1)>0)+(int_nu2x_cs(inxLL2)>0)+(int_nu2x_cs(inxLL3)>0)+(int_nu2x_cs(inxLL4)>0);
                else
                    inx_longerLL=zeros(l_seq1,1);
                end
                if (i<len-tot_l+1) && (tot_l<l) && half_nmerj<max_half_nmer
                    inxRRt=sum(seq_colp5n(inx_max_fit,i:i+nmer).*fivesi(ones(l_seq1,1),:),2)+1;
                    inxRR1=inxRRt+(1*(5^nmer));
                    inxRR2=inxRRt+(2*(5^nmer));
                    inxRR3=inxRRt+(3*(5^nmer));
                    inxRR4=inxRRt+(4*(5^nmer));
                    inx_longerRR=(int_nu2x_cs(inxRR1)>0)+(int_nu2x_cs(inxRR2)>0)+(int_nu2x_cs(inxRR3)>0)+(int_nu2x_cs(inxRR4)>0);
                else
                    inx_longerRR=zeros(l_seq1,1);
                end

                inx_nolong=find(inx_longerLL+inx_longerRR==0);
                int_chip_seq(inx_max_fit(inx_nolong))=max(int_chip_seq(inx_max_fit(inx_nolong)),int_nu2x_cs(inx(inx_max_fit(inx_nolong))));
                inx_long2=find(inx_longerLL|inx_longerRR>0);
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
                    int_chip_seq(inx_max_fit(inx_long2))=max(int_chip_seq(inx_max_fit(inx_long2)),intx);
                end
            end
        end
    end
end

end
