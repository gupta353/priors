% computation of marginal pdf of sigma from the joint Bernardo-pdf


function marginal_pdf_x=marginal_bernardo_pdf_sigma(x)

global GLOBAL_DATA
marginal_pdf_values_sigma=GLOBAL_DATA.marginal_pdf_values_sigma;
sigma_range=GLOBAL_DATA.sigma_range;

% compute marginal pdf at user-given x
diff=repmat(x,length(marginal_pdf_values_sigma),1)-sigma_range;
diff=abs(diff);
ind=find(diff==min(diff));

if length(ind)>1
    marginal_pdf_x=randomsample(marginal_pdf_values_sigma(ind),1);
else
    marginal_pdf_x=marginal_pdf_values_sigma(ind);
end

end