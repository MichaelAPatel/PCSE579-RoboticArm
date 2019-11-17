import numpy
import time
from scipy.integrate import odeint
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, show

# Define constants, note that moment of inertia will probably change.
mass = 0.05
radius = 0.6
g = 9.81
i = 2.0

# Method that approximates a second order ode.
def armmotion(vect, t):
    theta = vect[0]
    theta_dot = vect[1]
    tau = 5
    theta_ddot = (tau-(mass*g*radius*numpy.sin(theta))) / (i + pow(radius, 2))
    return [theta_dot, theta_ddot]

# Data from the two snapshots where angles are in degrees (Change as needed)
r1 = 5.0
theta1 = 15.8
r2 = 4.9
theta2 = 14.8

# Start the calculation clock.
start_time = time.time()

# TODO Insert trajectory calculations for the ball here. Once the trajectory is calculated we need to calculate the two
#  points it intersects the arms catchable region and the times at which those intersections occur. I presume we are
#  going to do the calculations in radians so just store the angle in the provided variable. Similarly, store the time
#  in the max_time variable.

angle = 'Insert Angle Here'
armangle = angle + (numpy.pi / 2)

max_time = 'Insert Time Here'

# Calculate arm motion by numerically solving the arm's equation of motion up to the max_time.
t = numpy.linspace(0.0, max_time, 101)
initialcondition = [0.0, 0.0]
solution = odeint(armmotion, initialcondition, t)

# Check when, if ever, the arm reaches the necessary angle.
for i in range(len(solution)):
    if solution[i][0] >= angle:
        arm_time = t[i]
        break

# Stop the calculation clock and determine if we are fast enough to catch the ball.
stop_time = time.time()
total_time = stop_time - start_time + arm_time

if total_time > max_time:
    print('The ball cannot be caught.')
else:
    print('The ball is caught at t=' + str(total_time) + 's.')


# This is some graphing stuff that we may want to use later.
# figure()
# plot(t, solution*180/numpy.pi)
# legend(('θ', 'θdot'))
# xlabel('TIME (sec)')
# ylabel('θ(deg) and θdot (deg/sec)')
# show()
