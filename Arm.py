import numpy
import time
import sys
import warnings
from scipy.integrate import odeint
from scipy.optimize import fsolve
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, show
warnings.filterwarnings('error')

# Define constants, note that moment of inertia will probably change.
mass = 0.03
radius = 0.6
g = 9.81
i = 0.0036


# Method to with equations that represent the balls trajectory and the catchable region. Used by FSolve to find
# intersection points.
def paths(z, x0, y0, vx, vy0):
    x = z[0]
    y = z[1]
    f = numpy.zeros(2)
    f[0] = pow(x, 2) + pow(y, 2) - pow(radius, 2)
    f[1] = y0 + vy0 * (x - x0) / vx - 4.905 * pow((x - x0) / vx, 2) - y
    return f


# Method that approximates a second order ode.
def armmotion(vect, t):
    theta = vect[0]
    theta_dot = vect[1]
    tau = -5
    theta_ddot = (tau - (mass * g * radius * numpy.sin(theta))) / (i + mass * pow(radius, 2))
    return [theta_dot, theta_ddot]


# Data from the two snapshots where angles are in degrees (Change as needed)
#Set 1
# r1 = 5.0
# theta1 = 225.13
# r2 = 3.8075
# theta2 = 225.32

#Set 2
r1 = 5.0
theta1 = 221.62
r2 = 3.8233
theta2 = 219.14

# Start the calculation clock.
start_time = time.time()

# Calculate X-Y positions and velocities based on the snapshot information.
x1 = r1 * numpy.cos(theta1 * numpy.pi / 180.0)
y1 = r1 * numpy.sin(theta1 * numpy.pi / 180.0)
x2 = r2 * numpy.cos(theta2 * numpy.pi / 180.0)
y2 = r2 * numpy.sin(theta2 * numpy.pi / 180.0)

vx = (x2 - x1) / 0.02
vy0 = (y2 - y1 - 0.001962) / 0.02

# Define initial guesses for the intersection points and some empty lists to hold intersection information.
intersectionGuesses = numpy.array([(radius, 0.0), (-1.0 * radius, 0)])
xint = []
tint = []
angleint = []

# Solve for intersection points, we solve twice using initial conditions on opposite sides of the catchable
# region that is defined by the arm's radius.
for i in intersectionGuesses:
    try:
        something = fsolve(paths, i, args=(x2, y2, vx, vy0))
        if something[0] not in xint:
            xint.append(something[0])
            angleint.append(numpy.arctan2(something[1], something[0]))
    except RuntimeWarning as e:
        next

# If there are no intersection points it is a ball and we can exit the program.
if not xint:
    stop_time = time.time()
    print('BALL!')
    print('\n')
    print('Initial X Position: ' + str(x2))
    print('Initial Y Position: ' + str(y2))
    print('Calcuated X Velocity: ' + str(vx))
    print('Calculated Initial Y Velocity: ' + str(vy0))
    print('Calculations completed in ' + str(stop_time - start_time) + 's.')
    sys.exit(0)

# Calculate the time for each intersection point based upon the x velocity which is constant.
for i in xint:
    tint.append((i - x2) / vx)

# Determine which angle we are going to use based on which one has the maximum time and then store that maximum time.
angle = angleint[tint.index(max(tint))]
max_time = max(tint)

# The arm angle is defined to be pointing straight down which is an angle pi/2 from the coordinates specified
# in the trajectory calculations. We have to convert the angle to arm angle and this is dependent upon which quadrant
# the original angle is in using symmetry arguments.
# First Quadrant
if 0 <= angle <= numpy.pi / 2:
    armangle = angle + (numpy.pi / 2)
# Second Quadrant
elif numpy.pi / 2 <= angle <= numpy.pi:
    armangle = 3 * numpy.pi / 2 - angle
# Third Quadrant
elif -1.0 * numpy.pi / 2 <= angle <= 0:
    armangle = angle + (numpy.pi / 2)
# Fourth Quadrant
elif -1.0 * numpy.pi <= angle <= -1.0 * numpy.pi / 2:
    armangle = numpy.abs(angle + (numpy.pi / 2))

# Calculate arm motion by numerically solving the arm's equation of motion up to the max_time.
t = numpy.linspace(0.0, max_time, 101)
initialcondition = [0.0, 0.0]
solution = odeint(armmotion, initialcondition, t)

# Check when, if ever, the arm reaches the necessary angle.
arm_time = pow(10, 17)
for i in range(len(solution)):
    if solution[i][0] >= armangle:
        arm_time = t[i]
        break

# Stop the calculation clock and determine if we are fast enough to catch the ball.
stop_time = time.time()
total_time = stop_time - start_time + arm_time

if total_time > max_time:
    print('STRIKE!')
else:
    print('OUT!')

print('\n')
print('Initial X Position: ' + str(x2))
print('Initial Y Position: ' + str(y2))
print('Calculated X Velocity: ' + str(vx))
print('Calculated Initial Y Velocity: ' + str(vy0))
print('Angles of intersection:')
for i in angleint:
    print('    ' + str(i * 180 / numpy.pi))
print('The program selected the angle ' + str(armangle * 180 / numpy.pi) + '. It must reach this angle in '
      + str(max_time) + ' to catch the ball. The arm was only able to reach ' + str(solution[len(solution)-1][0] * 180
      / numpy.pi) + '.')
print('Calculations completed in ' + str(stop_time - start_time) + 's.')

theta = []
for i in solution:
    theta.append(i[0]*180/numpy.pi)

figure()
plot(t, theta)
legend(('θ'))
xlabel('TIME (sec)')
ylabel('θ(deg)')
show()