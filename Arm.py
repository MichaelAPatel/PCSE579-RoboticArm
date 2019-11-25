import numpy
import time
import sys
import warnings
from scipy.integrate import odeint
from scipy.optimize import fsolve
from matplotlib.pyplot import figure, plot, show, grid, Circle, figtext
warnings.filterwarnings('error')

# ============================================Snapshot Information==================================================== #
# First Quadrant Out
r1 = 5.0
theta1 = 154.5
r2 = 4.8
theta2 = 153.86

# Second Quadrant Strike
# r1 = 5.0
# theta1 = 25.5
# r2 = 4.8
# theta2 = 26.14

# Third Quadrant Out
# r1 = 5.0
# theta1 = 25.5
# r2 = 4.8
# theta2 = 25.9

# Fourth Quadrant Out
# r1 = 5.0
# theta1 = 25.5
# r2 = 4.8
# theta2 = 25.7
# ==================================================================================================================== #

# Define constants, note that moment of inertia will probably change.
mass = 0.03
radius = 0.6
g = 9.81
i = 0.0036
tau = 5

annotateString = ''


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
    theta_ddot = (tau - (mass * g * (radius / 2) * numpy.sin(theta))) / (i + mass * pow(radius / 2, 2))
    return [theta_dot, theta_ddot]


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
    annotateString += 'BALL!'
    annotateString += '\nInitial X Position: ' + str(round(x2, 3))
    annotateString += '\nInitial Y Position: ' + str(round(y2, 3))
    annotateString += '\nCalculated X Velocity: ' + str(round(vx, 3))
    annotateString += '\nCalculated Initial Y Velocity: ' + str(round(vy0, 3))
    annotateString += '\nCalculations completed in ' + str(round(stop_time - start_time, 3)) + 's.'

    arm = Circle((0, 0), 0.6, color='b', fill=False)
    fig = figure(num=0.0, figsize=(6, 6))
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
    ax = fig.add_subplot(111, aspect='equal')
    grid(b=True, which='major', axis='both')
    ax.add_artist(arm)
    ax.set_xticks([-1.5, -1, -0.5, 0.5, 1, 1.5])
    ax.set_yticks([-1.5, -1, -0.5, 0.5, 1, 1.5])
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    figtext(0.0, 0.0, annotateString)
    show()
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
    armangle = -1.0 * (3 * numpy.pi / 2 - angle)
# Fourth Quadrant
elif -1.0 * numpy.pi / 2 <= angle <= 0:
    armangle = angle + (numpy.pi / 2)
# Third Quadrant
elif -1.0 * numpy.pi <= angle <= -1.0 * numpy.pi / 2:
    armangle = angle + (numpy.pi / 2)

# Calculate arm motion by numerically solving the arm's equation of motion up to the max_time.
t = numpy.linspace(0.0, max_time, 100)
initialcondition = [0.0, 0.0]
solution = odeint(armmotion, initialcondition, t)

# Check when, if ever, the arm reaches the necessary angle.
arm_time = pow(10, 17)
for i in range(len(solution)):
    if numpy.abs(solution[i][0]) >= numpy.abs(armangle):
        arm_time = t[i]
        break

# Stop the calculation clock and determine if we are fast enough to catch the ball.
stop_time = time.time()
total_time = stop_time - start_time + arm_time

if total_time > max_time:
    annotateString = 'STRIKE!'
else:
    annotateString = 'OUT!'

annotateString += '\n'
annotateString += '\nInitial X Position: ' + str(round(x2, 3))
annotateString += '\nInitial Y Position: ' + str(round(y2, 3))
annotateString += '\nCalculated X Velocity: ' + str(round(vx, 3))
annotateString += '\nCalculated Initial Y Velocity: ' + str(round(vy0, 3))
annotateString += '\nAngles of intersection:'
for i in angleint:
    annotateString += '\n    ' + str(round(i * 180 / numpy.pi, 3))
if total_time > max_time:
    annotateString += '\nThe program selected the angle ' + str(round(armangle * 180 / numpy.pi, 3)) + '. It must reach this angle in ' + str(round(max_time, 3)) + ' s. \n The arm was only able to reach ' + str(round(solution[len(solution) - 1][0] * 180 / numpy.pi, 3)) + ' in that time.'
else:
    annotateString += '\nThe program selected the angle ' + str(round(armangle * 180 / numpy.pi, 3)) + '. It must reach this angle in ' + str(round(max_time, 3)) + ' s. \n The arm was able to reach ' + str(round(armangle * 180 / numpy.pi, 3)) + ' in ' + str(round(arm_time,3)) + 's.'

annotateString += '\nCalculations completed in ' + str(round(stop_time - start_time, 3)) + 's.'

# Graph of Ball Trajectory
arm = Circle((0, 0), 0.6, color='b', fill=False)
x = numpy.linspace(-1, 1, 100)
y = y2 + vy0 * (x - x2) / vx - 4.905 * pow((x - x2) / vx, 2)
fig = figure(num=0.0, figsize=(6, 6))
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
ax = fig.add_subplot(111, aspect='equal')
grid(b=True, which='major', axis='both')
plot(x, y)
ax.add_artist(arm)
ax.set_xticks([-1.5, -1, -0.5, 0.5, 1, 1.5])
ax.set_yticks([-1.5, -1, -0.5, 0.5, 1, 1.5])
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
figtext(0.0, 0.0, annotateString)
show()
