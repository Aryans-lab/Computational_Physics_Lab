# Simulate the radioactive decay A → B → C for the decay constants λA = 0.0347
# and λB = 0.0198. Initial number of nuclei NA = 500, NB = 0, NC = 0. And plot overlay of Number of A, B and C versus time. (without using any libraries)
# Name: Aryan Bandyopadhyay, Roll number: 2411014

import matplotlib.pyplot as plt
from lcg import myrand

def radioactive_decay_simulation(NA0, NB0, NC0, lambda_A, lambda_B, t_max, dt):
    # Initials
    NA = [NA0]
    NB = [NB0]
    NC = [NC0]
    t = [0]

    n_steps = int(t_max / dt)

    # Pre-generate pool: NA0 * n_steps is safe since NA+NB <= NA0 always
    pool = [myrand() for _ in range(NA0 * n_steps)]
    idx = 0

    for step in range(n_steps):
        na = NA[-1]
        nb = NB[-1]

        # Slice random numbers for A and B nuclei this step
        ra = pool[idx : idx + na]
        idx += na
        rb = pool[idx : idx + nb]
        idx += nb

        # Count decays: nucleus decays if random number < lambda * dt
        decays_A = sum(1 for r in ra if r < lambda_A * dt)
        decays_B = sum(1 for r in rb if r < lambda_B * dt)

        NA.append(na - decays_A)
        NB.append(nb + decays_A - decays_B)
        NC.append(NC[-1] + decays_B)
        t.append(t[-1] + dt)

    return t, NA, NB, NC

t, NA, NB, NC = radioactive_decay_simulation(500, 0, 0, 0.0347, 0.0198, 400, 0.1)

plt.scatter(t, NA, s=2, label='NA')
plt.scatter(t, NB, s=2, label='NB')
plt.scatter(t, NC, s=2, label='NC')
plt.xlabel('Time')
plt.ylabel('Number of nuclei')
plt.title('Radioactive decay A → B → C')
plt.legend()
plt.tight_layout()
plt.show()

# End of Code