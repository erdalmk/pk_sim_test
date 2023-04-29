# pk_sim_test
Simulate reasonable pharmacokinetic model using SimBiology and fit noisy data using IPOPT through CasADi.

Simbiology produces nonuniform timestamps, so we use exact discretizations which is provided by two symbolicly created functions. These functions will be created in the first run and if they are in the running directory, they will be directly used. Hence the first run will be slower.
