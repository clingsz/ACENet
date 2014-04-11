function [x] = reshapeTwins(X,familyId)
% reshapeTwins takes an twin based arrays of measurements and
% family ids and returns family based arrays of measurements
%
% [x,z] = reshapeTwins(X,familyId) 
%    X        - an array of size Phenos x Twins contains phenotype
%               measurements for each twin
%    familyId - a vector of family ids for each Twin, two twins from the
%               same family have the same familyId
% Returns
%    x        - an array of size 2 x Families x Phenos contains phenotype
%               measurements for each family (2 twins)


[cs,ps] = countDistinct(familyId);
ps = ps(find(cs==2));

lst = findIn(familyId,ps);
[~,ii] = sort(familyId(lst));
lst = lst(ii);

x = zeros(2,size(X,2)/2,size(X,1));

res2 = @(v) reshape(v,[2 length(v)/2]);
for i=1:size(X,1)
    x(:,:,i) = res2(X(i,:));
end



