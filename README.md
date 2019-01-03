# SubCom v1.0 
#(SUBsidence due to COMpaction of Aquifers version 1.0)

1. General information.

SubCom v1.0 program was implemented in the Scilab environment (version 5.5.2). It can be used for 
determining model parameters of compaction layers (example in oil, gas or groundwater extraction). 
Program solved the inverse problem and calculated: compaction coefficient Cm of reservoir rocks mass 
and parameter tgB, which characterizes the properties of the overburden.

2. Ho to use?
    
    a. Run the 'SubCom_v1_0.sce' script (recommended Scilab v5.5.2).The GUI will be displayed.
	
    b. In main windows program set the path to Aquifer and Profile information. Set also the initial 
value of parameters and number of iterations.

    c. Check the corecness of Aquifer and Profile data by 'Draw data' buttom.
	
    d. Click 'Calculations' buttom. Model parameters for the compaction reservoir based on the 
least square method will be calculated and the result will be presented systematically on the screen.

    e. The results will be written to the 'Results.txt' file in main folder.
	
    f. The '1-Iteration' buttom calculated teretical value of subsidence based on curent value of Cm and tgB.
	
The result will be ploted with empirical profile.

Contact and Information: Wojciech Witkowski (e-mail: wwitkow@agh.edu.pl) for questions regadring 
the Scialab codes.
