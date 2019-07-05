% computation of marginal pdf of kh from the joint Bernardo-pdf


function marginal_pdf_x=marginal_bernardo_pdf_kh(x)

global GLOBAL_DATA
marginal_pdf_values_kh=GLOBAL_DATA.marginal_pdf_values_sigma;
kh_range=GLOBAL_DATA.kh_range;

% compute marginal pdf at user-given x
diff=repmat(x,length(marginal_pdf_values_kh),1)-kh_range;
diff=abs(diff);
ind=find(diff==min(diff));

if length(ind)>1
    marginal_pdf_x=randomsample(marginal_pdf_values_kh(ind),1);
else
    marginal_pdf_x=marginal_pdf_values_kh(ind);
end

end