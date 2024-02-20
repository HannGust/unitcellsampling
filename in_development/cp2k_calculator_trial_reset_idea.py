#

# Contains the code for a trial start/reset portion of grid_gen_script.py
# to see if it could be a work around for some probelms with the cp2k shell,
# and the MPI/OMP problems on the cluster.

# It worked, but didn't solve the issue as it turned out to be related to something else.

##### Perhaps test: Try - except statement here to see if it catches failure to start error
##### debbuging test:
cp2k_calc_not_started = True
cp2k_calc_start_attempts = 1
calc_start_limit = 10

#### Lets try to overwrite the cp2kshell expect method:

while cp2k_calc_start_attempts <= calc_start_limit  and cp2k_calc_not_started:
    try:
        cp2k_calc_from_inp = cp2k_calculator_from_input(parsed_cp2k_input, 
                                                    **cp2k_kwargs)
    except Exception as exc:
        print("Starting cp2k calc raised exception: ", str(exc))
        print("Type: ", type(exc))
        print("Args: ", exc.args)
        print("Dict: ", exc.__dict__)
        print("Cause: ", exc.__cause__)
        print("Context: ", exc.__cause__)
        print("Start attempt: ", str(cp2k_calc_start_attempts))
        cp2k_calc_start_attempts += 1
    else:
        print("Cp2k calculator successfully started. Attempt: ", str(cp2k_calc_start_attempts))
        cp2k_calc_not_started = False

if (cp2k_calc_not_started) and (cp2k_calc_start_attempts > calc_start_limit):
    print("Could not start cp2k calculator. Aborting.")
    raise Exception("Could not start cp2k caluclator after multiple attempts. Aborting.")
 
####

