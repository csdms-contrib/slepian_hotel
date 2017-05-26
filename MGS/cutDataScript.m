% Reduce all data once you have the indices selected

Bdyn{1}=Bdyn{1}(indices);
Bdyn{2}=Bdyn{2}(indices);
Bdyn{3}=Bdyn{3}(indices);

Brms{1}=Brms{1}(indices);
Brms{2}=Brms{2}(indices);
Brms{3}=Brms{3}(indices);

Bstat{1}=Bstat{1}(indices);
Bstat{2}=Bstat{2}(indices);
Bstat{3}=Bstat{3}(indices);

cola=cola(indices);
lon=lon(indices);

data{1}=data{1}(indices);
data{2}=data{2}(indices);
data{3}=data{3}(indices);

r=r(indices);

sam=sam(indices);
sao=sao(indices);
sap=sap(indices);

date=date(indices,:);