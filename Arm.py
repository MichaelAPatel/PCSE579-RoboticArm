import numpy
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, show
from scipy.integrate import odeint


mass = 0.05
radius = 0.6
g = 9.81
i = 2.0

def armmotion(vect, t):
    theta = vect[0]
    theta_dot = vect[1]
    tau = 5
    theta_ddot = (tau-(mass*g*radius*numpy.sin(theta))) / (i + pow(radius, 2))
    return [theta_dot, theta_ddot]


t = numpy.linspace(0.0,1.0,101)
initialcondition = [0.0, 0.0]
solution = odeint(armmotion, initialcondition, t)

figure()
plot(t, solution*180/numpy.pi)
legend(('θ', 'θdot'))
xlabel('TIME (sec)')
ylabel('θ(deg) and θdot (deg/sec)')
show()
