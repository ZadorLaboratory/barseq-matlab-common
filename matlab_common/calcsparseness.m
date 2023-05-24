function s=calcsparseness(x)
%calculate the sparseness of x. x is a n x m matrix with n rows indicating
%samples and m columns indicating bins. sparseness is calculated as
%1-(sum(x0^2)/N)/(sum(x0)/N)^2
s=1-((sum(x,2)/size(x,2)).^2)./(sum(x.^2,2)/size(x,2));
s(isnan(s))=1;
end
