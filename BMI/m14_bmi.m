% docker run -ti --rm -v $PWD:/opt/MARRMoT --entrypoint bash sverhoeven/marrmot
% cd /opt/MARRMoT
% octave

addpath(genpath('.'))
model = marrmotBMI_oct()
model.initialize('BMI/Config/BMI_testcase_m14_BuffaloRiver_TN_USA.mat')
model.update()
model.get_value('flux_out_Q')
