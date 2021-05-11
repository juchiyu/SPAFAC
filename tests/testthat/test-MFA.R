library(mpMFA)

FactoMineR::

data('wines2007')
demo.mfa.2007 <- mpMFA(wines2007$data, wines2007$table, graphs = FALSE)

X = wines2007$data
column.design = sub("\\..*", "", colnames(wines2007$table))
components = 2L
sparseOption = "subtable";
center = TRUE; mfa.scale = TRUE;
tol = .Machine$double.eps;
init = "svd"; initLeft = NULL; initRight = NULL; seed = NULL;
rdsLeft = rep(0.7*sqrt(nrow(wines2007$data)), components); rdsRight = rep(1, components);
grpLeft = NULL;
orthogonality = "both";
OrthSpaceLeft = NULL; OrthSpaceRight = NULL;
projPriority = "orth";
projPriorityLeft = projPriority;
projPriorityRight = projPriority;
itermaxALS = 1000; itermaxPOCS = 1000;
epsALS = 1e-10; epsPOCS = 1e-10
