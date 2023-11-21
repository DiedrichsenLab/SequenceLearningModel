import numpy as np

def do_sim(M, T, dT=2, maxTime=200000,plan_func='plan_exp',trigger_func='trig_indep'):

    T.reinit()
    T.stimTime = np.full(T.numPress, np.nan)
    T.futureTime = np.full((T.numPress,M.wFuture), np.nan) # time that triggering horizon is reached for each element

    # Initialize variables
    X = np.zeros((M.numOptions, maxTime // dT, T.numPress))  # Hidden state
    S = np.zeros((M.numOptions, maxTime // dT, T.numPress))  # Stimulus present
    B = np.ones(maxTime // dT) * M.Bound  # Decision Bound - constant value
    t = np.arange(1, maxTime // dT + 1) * dT - dT  # Time in ms


    i = 0  # Index of simulation
    nDecision = 1  # Current decision to make
    numPresses = 0  # Number of presses produced
    isPressing = 0  # Is the motor system currently occupied?

    # Set up parameters
    dec = np.arange(1, T.numPress + 1)  # Number of decision
    A = np.eye(M.numOptions) * (M.Aintegrate - M.Ainhibit) + np.ones((M.numOptions, M.numOptions)) * M.Ainhibit

    # Start time-by-time simulation
    while numPresses < T.numPress and i < maxTime // dT:

        # horizon for triggering rule
        horizon = np.arange(nDecision,T.numPress)[:M.wFuture]

        # Figure out when stimuli appear
        visible = np.arange(min(T.numPress, numPresses + T.window))
        T.stimTime[visible[np.isnan(T.stimTime[visible])]] = t[i]

        # Update the stimulus: Fixed stimulus time
        indx = np.where(t[i] > (T.stimTime + M.dT_visual))[0]  # Index of which stimuli are present T.
        if len(indx) > 0:
            for j in indx:
                S[T.stimulus[j] - 1, i, j] = 1

        # Determine the distribution of rates planning
        rate = getattr(M, plan_func)(dec,nDecision)

        # Update the horse race
        eps = np.random.randn(M.numOptions, 1, T.numPress) * M.SigEps  # Noise
        for j in range(T.numPress):
            X[:, i + 1, j] = np.dot(A, X[:, i, j]) + rate[j] * S[:, i, j] + dT * eps[:, 0, j]

        # Determine if we issue a decision
        if nDecision <= T.numPress and not isPressing and np.any(X[:, i + 1, nDecision - 1] > B[i + 1]):
            
            # check if all future decisions are above bound
            horizoncheck, element_done = getattr(M, trigger_func)(X[:,i+1,horizon], rate[horizon])
            mask = np.isnan(T.futureTime[nDecision-1, :]) & (element_done == 1)
            T.futureTime[nDecision-1,mask] = t[i+1]

            if horizoncheck and (nDecision > 1 or (T.RT is None) or (T.RT is not None and t[i] > T.RT)):
                T.response.append(np.argmax(X[:, i + 1, nDecision - 1])+1)
                T.decisionTime.append(t[i + 1])  # Decision made at this time
                T.pressTime.append(T.decisionTime[nDecision - 1] + M.dT_motor)  # Press time delayed by motor delay

                isPressing = 1  # Motor system engaged
                nDecision += 1

        # Update the motor system: Checking if current movement is done
        if isPressing and t[i + 1] >= T.pressTime[nDecision - 2]:
            isPressing = 0
            numPresses += 1

        i += 1

    SIM = {'X': X[:, :i, :], 'S': S[:, :i, :], 'B': B[:i], 't': t[:i]}
    return T,SIM