import numpy as np
x=np.load("./solo_out/softmax_scores.npy")
np.savetxt('solo.txt',x)
