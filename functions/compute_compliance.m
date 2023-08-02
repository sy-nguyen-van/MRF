function [c,grad_c] = compute_compliance()

global FE OPT
% compute the compliance
c = (FE.P)'*FE.U;
Ue = permute(repmat( FE.U(FE.edofMat).', [1,1,FE.n_edof]), [1,3,2]);
Ue_trans = permute(Ue, [2,1,3]);%
grad_c = reshape(sum(sum(-Ue_trans.*FE.dKe.*Ue, 1),2),[1,FE.n_elem]);
grad_c = OPT.H*((OPT.dv*1.875).*grad_c(:))./OPT.Hs./(OPT.dv*1.875); % 99 line codes: EQ. 5   
%------------------------------------------------------------
OPT.compliance = c;
OPT.grad_compliance = grad_c;
end
