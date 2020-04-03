#!/usr/bin/python


#######################################################
# Python Script to plot rms values (for theta, phi, y)
# for different minimisers and algorithms
#  John Williamson
#  21/1/2019
#######################################################


import csv



def extract_var(var, csv_reader):

    var_results = []


    print('EXTRACT_VAR !!! \n')
    for row in csv_reader:

        if(row['Optimisation Variable']==var):
           var_results.append(row)


    print('EXIT EXTRACT_VAR !!! \n')
    return var_results


with open('rms_csv/algorithm_results_toy.csv') as csv_file:
    csv_Dict = csv.DictReader(csv_file, delimiter=',')

    csv_reader = list(csv_Dict)
    #    line_count = 0
    



    import matplotlib.pyplot as plt
    
    #set size of image
    plt.rcParams["figure.figsize"] = (20,5)

    
    cm = plt.get_cmap('gist_rainbow')

    fig = plt.figure("Comparison of Minimisers and algorithms")


    ## set-up loop to go through different optimisation variables (theta,phi,y) and plot results

    opt_vars = ['Theta','Phi','Y']


    # set-up list of markers to iterate through for different minimisers
    markers = ["o","v","s"]


    var_int = 0 # iterator to track opt_var


    for var in opt_vars:

        print('\n new var: ',var, '\n')
        
        var_int +=1



        var_results = []


        var_results = extract_var(var,csv_reader)
        

        RMS = [rms['Rms (real xyz)'] for rms in var_results]
        RMS = [float(i) for i in RMS] # line converts from string to float

        Minimiser = [mini['Minimiser'] for mini in var_results]

        Algorithm = [alg['Algorithm'] for alg in var_results]
        


        

        Min_set = set(Minimiser)
        
        all_len = len(Minimiser)
        print('type(max(RMS)) = ', type(max(RMS)))
        

    
        var_plot = fig.add_subplot(3,1,var_int)

        #        print(max(RMS))
        maxRMS = max(RMS)
        minRMS = min(RMS)
        rangeRMS = maxRMS - minRMS

        var_plot.set_ylim(minRMS-rangeRMS,maxRMS+rangeRMS)
        

        var_plot.set_ylabel(var + ' RMS')
        var_plot.set_xlabel('Minimiser')
        #        var_plot.set_ylim([0,1])

        h = 0 # iterator for different minimisers

        for i in Min_set:

        
            print('for ', i)
        
            k = 0 # iterator for all entries in table
            print('k = ', k)


            col_len = len([a for a in Minimiser if a == i])

            print('col_len = ', col_len)

            
            l = 0 # iterator for differnt algorithms of one minimiser

            for j in Minimiser:
            
        
                if j == i:
                    print(k)
                    #point = plt.plot(j,RMS[k],markers[h], label = Algorithm[k])
                    point = var_plot.plot(j,RMS[k],markers[h], label = Algorithm[k])
                    print("l = ", l)
                    #                point[0].set_color(cm(l//3*3.0/col_len)) 
                    point[0].set_color(cm(l*1.0/col_len)) 
                    # point[0].set_color(cm(k*1.0/all_len)) 
                    l += 1
        
                k += 1
            print('\n')                        
            h +=1
        
        
        if var == 'Theta':
            var_plot.legend(loc = 'lower left', bbox_to_anchor= (0.0, 1.01), ncol=2)
        


#    var_plot.legend()
    

    
    plt.show()
