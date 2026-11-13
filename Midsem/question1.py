# Diffusion of particles between two halves of a box
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

import matplotlib.pyplot as plt
from mylib import myrand

def diffusion(N, steps):
    # Track particle locations: 0 for Left, 1 for Right
    particles = [0] * N
    N_left = N
    N_right = 0

    time_steps = [0]
    N_left_history = [N_left]
    N_right_history = [N_right]

    for t in range(1, steps + 1):
        # Pick one particle at random from all N particles
        i = int(myrand() * N)
        if particles[i] == 0:
            particles[i] = 1
            N_left -= 1
            N_right += 1
        else:
            particles[i] = 0
            N_left += 1
            N_right -= 1

        time_steps.append(t)
        N_left_history.append(N_left)
        N_right_history.append(N_right)

    return time_steps, N_left_history, N_right_history

if __name__ == "__main__":
    N = 5000  # Total particles
    steps = 30000  # Should be Sufficient time (>> N = 5000)
    time_steps, N_left, N_right = diffusion(N, steps)

    # Plotting number of particles on left and right sides vs time
    plt.plot(time_steps, N_left, label='Left side', color='blue')
    plt.plot(time_steps, N_right, label='Right side', color='orange')
    plt.axhline(N / 2, color='black', linestyle='--', label=f'Equilibrium (N/2 = {N // 2})')
    plt.scatter([time_steps[-1]], [N_left[-1]], color='blue')
    plt.scatter([time_steps[-1]], [N_right[-1]], color='orange')
    plt.annotate(f'Left: {N_left[-1]}', (time_steps[-1], N_left[-1]), textcoords="offset points", xytext=(-70, 10))
    plt.annotate(f'Right: {N_right[-1]}', (time_steps[-1], N_right[-1]), textcoords="offset points", xytext=(-70, -15))
    plt.xlabel('Time steps')
    plt.ylabel('Number of particles')
    plt.title('Diffusion of particles between two halves of a box')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Equilibrium state determination from graph
    print(f"Total particles N = {N}")
    print(f"At final time step t = {steps}:")
    print(f"Number of particles on left side  = {N_left[-1]}")
    print(f"Number of particles on right side = {N_right[-1]}")
    print(f"Both sides reach equilibrium around N/2 = {N // 2}")

# End of code

##########################################################
# Output:
# Total particles N = 5000
# At final time step t = 30000:
# Number of particles on left side  = 2468
# Number of particles on right side = 2532
# Both sides reach equilibrium around N/2 = 2500
###########################################################