% PROGRAM TO SCORE SEQUENCES USING MINSEQS - v2
function int_chip_seq=...
    chip_seq_count_v2(max_gap,min_half_nmer,max_half_nmer,seqsx_inx_c,l,...
    pow_factor,all_len_nox,int_nu2x_c,bed_file,min_nmer,max_nmer)

if max_half_nmer==-1
    m_size=5^max_nmer;
else
    max_gap=min(max_gap,l-2*min_half_nmer);
    m_size=(5^(2*max_half_nmer))*(max_gap+1);
end

int_nu2x_cx=int_nu2x_c.*(pow_factor.^all_len_nox);
seqsx_inxu_rev_c=get_rev_inx_v2(seqsx_inx_c,max_half_nmer,l);
[seq_uq,ia,~]=unique([seqsx_inx_c seqsx_inxu_rev_c]);
int_nu2x_cx2=[int_nu2x_cx int_nu2x_cx];
int_nu2x_cs=sparse(seq_uq,1,int_nu2x_cx2(ia),m_size,1);
tot_chips=1;
int_chip_seq=cell(tot_chips,1);

for jj=1:tot_chips
    ip_file=bed_file;
    C=textscan_mod_v1(ip_file,'%s','\t');
    seqs_chip=C{1};
    cont_seqs=zeros(length(seqs_chip),1);
    for j=1:length(seqs_chip)
        if length(seqs_chip{j})<1000
            cont_seqs(j)=1;
        end
    end
    ix=find(cont_seqs==1);
    seqs_chip=seqs_chip(ix);
    [seq_colp5n_chip,~,l_seqs_chip]=get_seq_col_v3(seqs_chip,[],[]);
    max_len=0;
    for i=1:l_seqs_chip
        max_len=max(max_len,length(seqs_chip{i}));
    end
    seq_colp5n_chipx=zeros(l_seqs_chip,max_len);
    for i=1:l_seqs_chip
        seq_colp5n_chipx(i,1:length(seqs_chip{i}))=seq_colp5n_chip{i};
    end

    if max_half_nmer==-1
        int_chip_seq{jj}=get_seq_scores_ng_v2(seq_colp5n_chipx,l,...
            min_nmer,max_nmer,int_nu2x_cs);
    else
        int_chip_seq{jj}=get_seq_scores_gp_v2(seq_colp5n_chipx,l,...
            min_half_nmer,max_half_nmer,max_gap,int_nu2x_cs);
    end
end


end
