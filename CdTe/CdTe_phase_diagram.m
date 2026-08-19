

T_Cdrich = [593.9 1123 1198 1245 1283 1317 1350 1365];
X_Cdrich = [0 5.021 15.001 25 34.986 42.493 47.456 50.089];

T_Terich = [1365 1347 1315 1250 1198 1148 1088 977 723];
X_Terich = [50.089 51.996 54.994 60 65.235 72.25 80 89.989 100];


X_poly_Cdrich = polyfit(T_Cdrich,X_Cdrich,numel(X_Cdrich));  % give x for an input T
X_poly_Terich = polyfit(T_Terich,X_Terich,numel(X_Terich));

fit_Cdrich = polyval(X_poly_Cdrich, T_Cdrich);
fit_Terich = polyval(X_poly_Terich, T_Terich);

figure
clf
hold on
plot(T_Cdrich,X_Cdrich,'bo')
plot(T_Terich,X_Terich,'ro')
plot(T_Cdrich,fit_Cdrich,'b-')
plot(T_Terich,fit_Terich,'r-')



