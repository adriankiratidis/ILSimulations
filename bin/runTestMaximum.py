#### Shell script to run runSingleSphere.csh
#!/usr/bin/python
import os

min_try=20 #Suggested to set to largest separation of previous success
max_try=40 #Suggested to set to smallest separation of previous failure

C_attempted=[]

#Start with some numbers so our algorithm can take the difference.
succesfull_attempts=[0] #Zero should always work
failed_attempts=[1000] #very large separations should fail

tolerance = 1 #The difference between the largest succesful and smallest failed attempt
#
while True:
    
    #Set a value of separation
    current_min_attempt = max(min_try, max(succesfull_attempts))
    current_max_attempt = min(max_try, min(failed_attempts))
    C = (current_min_attempt + current_max_attempt)/2

    if(C in C_attempted):
        print("Tried all possible C's to extremums....exiting....")
        print("C_attempted = ", C_attempted)
        exit()
    
    #Run at that separation
    C_attempted.append(C)
    print("Running C = ", C)
    os.system("bash bin/runSimulation15.sh {0}".format(str(C)))

    #Did I run succesfully?
    output = os.popen("ls /home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_PositiveMinusSpheres-0.005-" + str(C) + "--0.003125-hs1.0TESTING_MAXIMA_FINDER | grep 'testing15-a2-potential-per-unit-area.txt'").read()
    #print("output = ", output)

    if(output.strip() == ''):
        print("Ran UNSUCCESSFULLY for C = ", C)
        failed_attempts.append(C)
    else:
        print("Ran SUCCESSFULLY for C = ", C)
        succesfull_attempts.append(C)


    #Check to make sure there is no overlap between succesful and unsuccesful tries.
    if(max(succesfull_attempts) > min(failed_attempts)):
        print("max(succesfull_attempts) > min(failed_attempts)")
        print("This case is not dealt with by the script")
        exit()

    #Check if we've found maxima within tolerance, if so, exit
    if(min(failed_attempts) - max(succesfull_attempts) == tolerance):
        print("completed finding maxima run")
        print("maxima = ", max(successfull_attempts))
        exit()
