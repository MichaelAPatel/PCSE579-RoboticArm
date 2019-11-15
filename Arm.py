import numpy
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, show
#from scipy.integrate import odeint
from math import sin, cos, pi
from busect import insort
mass = 0.05
radius = 0.6
g = 9.81
i = 2.0
armlink = .3
forearmlink = .3
catchrad = .075  #this should be the ball radius plus the hand radius or how far away the center of the hand needs to be from the center of the ball to catch it(.075 diameter of baseball placeholder)

def dvSholder (T,arm, dt):
    return 1
def dvElbow (T,arm, dt):
    return 1
def dTh(angV, dt):
    return angV * dt
class ArmStruct:
    def __init__(self, lastStruct, sholderT, elbowT, timestep = .01):
        if lastStruct:
            self.time = lastStruct.time + timestep
            self.sholderTh = lastStruct.sholderTh + dTh(lastStruct.sholderV, timestep)
            self.elbowTh   = lastStruct.sholderTh + dTh(lastStruct.elbowV,   timestep)
            self.sholderV = lastStruct.sholderV + dvSholder(sholderT, lastStruct, timestep)
            self.elbowV = lastStruct.sholderV + dvSholder(sholderT, lastStruct, timestep)
            self.sholderT = sholderT
            self.elbowT = elbowT
            self.handX = armlink * cos(self.sholderTh) + forearmlink * cos(self.sholderTh + self.elbowTh)
            self.handY = armlink * sin(self.sholderTh) + forearmlink * sin(self.sholderTh + self.elbowTh)
        else:
            self.time = 0
            self.sholderTh = - pi/2
            self.elbowTh = 0
            self.sholderV = 0
            self.elbowV = 0
            self.sholderT = 0
            self.elbowT = 0
            self.handX = armlink * cos(self.sholderTh) + forearmlink * cos(self.sholderTh + self.elbowTh)
            self.handY = armlink * sin(self.sholderTh) + forearmlink * sin(self.sholderTh + self.elbowTh)


def genSuccessors(arm):
    """
    this will create a list of arm structs that can be acheived from the current state at this time step. I will create
    successors that will have at leat one of the acuators at a torque of 5 or -5  and the other with 11 other torques making up to
    29 possibilities in total.
    :param arm: this is an arm stuct of an arm that is in the inital position
    :return: a list of arm structs of possible successors after a timestep
    """
    succlist = []
    for eT in range(-5,6,1):
        for sT in(-5,5):
            succlist.append(ArmStruct(arm, sT,eT)
    for sT in range(-4,5,1):
        for eT in (-5,5):
            succlist.append(ArmStruct(arm, sT, eT)
    return list(filter(lambda x: x.sholderV < 60*pi and x.elbowV < 60*pi, succlist)

def isGoal(handposition,goalstates):
    """

    :param handposition: this is a tuple of floats the (x,y, time) that describe the locaton of the hand in 2d*time space
    :param goalstates: this is a tuple of goal states of the same syntax as hand positions that are in the reach of the ball
    :return: true if the hand position is already able to catch the ball and false otherwise
    """
    for goal in goalstates:
        if handposition[2] < goal[2]:
            if ((handposition[0] - goal[0])**2 + (handposition[1] - goal[1])**2)**0.5 < catchrad:
                return True
    return False
def heu(handposition,goalstates):
    """

    :param handposition: this is a tuple of floats the (x,y, time) that describe the locaton of the hand in 2d*time space
    :param goalstates: this is a tuple of goal states of the same syntax as hand positions that are in the reach of the ball
    :return: an optomistic heuristic for the getting to the most optomistic goal state
    """
    best = float("inf")
    for goal in goalstates:
        if handposition[2] < goal[2]:
            best = min(best,(((handposition[0] - goal[0]) ** 2 + (handposition[1] - goal[1]) ** 2) ** 0.5) /(goal[2] - handposition[2]))

    return best


def bisectInsert(item, sortedList, key = lambda x: x):

    """
    this does in place insteriton of the item into the list to stay sorted
    :param item: item to insert into list
    :param inPlaceList: the list that will be inserted to
    :param key: this decides what to consider while sorting the list
    :return: a sorted list with the new item inserted
    """
    low = 0
    high = len(list)
    while True:
        if low == high:
            return sortedList[:low]+[item]+ sortedList[high:]
        ndx = int((high - low)/2)+ low
        if compare(item)> compare(sortedList[ndx]):
            low = ndx +1
            continue
        high = ndx


def solutionSearch(initialArm, goalstates):
    """
    this will do a greedy search for a catch condition of the ball using a graph search
    :param initialArm: this i the initial state of the arm
    :param goalstates:
    :return:
    """
    fringe =[(0,initialArm,[])]
    while fringe:
        _, thisState, thispath = fringe.pop()
        if isGoal((thisState.handX, thisState.handY, thisState.time), goalstates):
            return thisState, thispath
        succs = genSuccessors(thisState)
        succs = [(heu(x), x, thispath+[x]) for x in succs]
        succs = list(filter(labda x: x[0] != float("inf"), succs))
        for item in succs:
            fringe = bisectInsert(item, fringe, key= lambda x: x[0])



def armmotion(vect, t):
    theta = vect[0]
    theta_dot = vect[1]
    tau = 5
    theta_ddot = (tau-(mass*g*radius*numpy.sin(theta))) / (i + pow(radius, 2))
    return [theta_dot, theta_ddot]


t = numpy.linspace(0.0, 1.0, 101)
initialcondition = [0.0, 0.0]
solution = odeint(armmotion, initialcondition, t)

figure()
plot(t, solution*180/numpy.pi)
legend(('θ', 'θdot'))
xlabel('TIME (sec)')
ylabel('θ(deg) and θdot (deg/sec)')
show()
