% FUNCTION seq_col & seqs_mat MATRIX FROM SEQUENCES v3 (Sequences of different lengths)
function [seqs_mat,seq_col,l_seqs]=get_seq_col_v3(seqs,Llib,Rlib)
aa2=[1 0 2 0 0 0 3 0 0 0 0 0 0 5 0 0 0 0 0 4 0 0 0 0 0 0];
aa2(56)=5; %FOR x
nucl_to_mat=[eye(4) ; ones(1,4)]; % FOR N
if iscell(seqs)
    l_seqs=length(seqs);
    l=length(seqs{1});
    seqs_mat=cell(l_seqs,1);
    seq_col=cell(l_seqs,1);
    for j=1:l_seqs
        seqs_mat{j}=[aa2(Llib-64) aa2(seqs{j}-64) aa2(Rlib-64)];
        seq_col{j}=sparse(reshape(nucl_to_mat(seqs_mat{j},:)',length(seqs_mat{j})*4,1)');
    end

else
    [l_seqs,l]=size(seqs);
    seqs_mat=[repmat(aa2(Llib-64),[l_seqs 1]) aa2(seqs-64) repmat(aa2(Rlib-64),[l_seqs 1])];
    seq_col=reshape(nucl_to_mat(seqs_mat',:)',(l+length(Llib)+length(Rlib))*4,l_seqs)';
end

end
