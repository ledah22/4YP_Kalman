# 4YP_Kalman
Welcome to the repository of my 4YP project "Control of Transcranial Magnetic Stimulators: Kalman Filter based approach to timed stimulation".  
Folder *ForExaminers* is to be used for reproducing results presented in the Results section of the report.  
Code writing was approached in a modular manner, so systems for artefact removal, phase tracking and the combination of the two can be tested separately.

To test the final product, please run *Kalman_Modular.m* script to set the parameters. Then run the corresponding Simulink model in *Kalman_Modular_Simple.slx*. 
In the Matlab script, it is possible to set the type of test data (synthesised or accelerometer tremor data), but the accelerometer data is only available upon request.  
Different testing setups can also be experimented with: different frequencies, length of simulation, number of samples used to tune parameters etc.

To explore individual modules, use *Kalman_EMG.m* and *Kalman_Simple.slx* for the artefact removal algorithm.  
Run *Kalman_Phase.m* and *Kalman_Phase_Simple.slx* to test the phase tracker.

Finally, there are functions available to plot the results and compare them to the Hilbert transform and input signal.  
* *make_plot.m* was primarily designed to visualise the results of the artefact removal algorithm 
* *phaseplot.m* displays the outputs of tracking phase  
* *combinedplot.m* visualises the outputs of the combined system and can also provide error calculations 
* *tremorplot.m* was tailored to clearly display the outcomes of analysing real accelerometer data 
