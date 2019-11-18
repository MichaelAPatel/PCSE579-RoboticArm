import numpy
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, show
#from scipy.integrate import odeint
from math import sin, cos, pi
mass = 0.05
radius = 0.6
g = 9.81
i = 2.0
armlink = .3
forearmlink = .3
catchrad = .075  #this should be the ball radius plus the hand radius or how far away the center of the hand needs to be from the center of the ball to catch it(.075 diameter of baseball placeholder)


def findGoalsFromStart(point1, point2, dt = .002):
    """
    this function will return various points where the ball will be located inside of the reach of the arm. These points
    should allow for overlapping of the ball from one point to the next but not so that they are overlapping a lot.
    :param point1: this will be a tuple of the form (float, float) of the point at time t= -dt in radial coordinates with
    the first float being the distance to the origin and the second float being the angle between the location of the
    ball and 0 (with 0 being to the right)
    :param point2: same as point 1 just at time t= 0
    :return: this should return a list of tuples of the form (float, float, float) where the first float is the x
    coordinate of a location along the ball path, the sencond is the y coordinate and the last is the time coordinate.
    If there are no goal states return an Empty list
    """
    goals = []
    #########PUT CODE HERE###########

    return goals


def dvShoulderAndElbow (T, arm, dt):

    """
    called when initalizing an Arm Struct used to find the new anglular velocity of the shoulder acuator
    :param T: a tuple of (float, float)  that represents the torque of the shoulder for [0] and elbow at [1]
    :param arm: this is a ArmStruct that represents the last position of the arm
    :param dt: this is a float that represents the time step that you have taken since your last known theta
    :return: a tuple of (float, float)  that represents the angualr velocity of the shoulder at the next time step for [0]
    and elbow at [1]
    """
    return 1., 1.


def dTh(angV, dt):
    """
    called when initalizing an Arm Struct used to find the new angle
    :param angV: is a float that represents angular velocity of your system
    :param dt: this is a float that represents the time step that you have taken since your last known theta
    :return: the change in theta since your last timestep
    """
    return angV * dt


class ArmStruct:
    def __init__(self, lastStruct, shoulderT, elbowT, timestep = .01):
        if lastStruct:
            self.time = lastStruct.time + timestep
            self.shoulderTh = lastStruct.shoulderTh + dTh(lastStruct.shoulderV, timestep)
            self.elbowTh =    lastStruct.elbowTh +    dTh(lastStruct.elbowV,    timestep)
            sdV, edV = dvShoulderAndElbow((shoulderT, elbowT), arm, dt)
            self.shoulderV = lastStruct.shoulderV + sdV
            self.elbowV =    lastStruct.elbowV    + edV
            self.shoulderT = shoulderT
            self.elbowT = elbowT
            self.handX = armlink * cos(self.shoulderTh) + forearmlink * cos(self.shoulderTh + self.elbowTh)
            self.handY = armlink * sin(self.shoulderTh) + forearmlink * sin(self.shoulderTh + self.elbowTh)
        else:
            self.time = 0
            self.shoulderTh = - pi/2
            self.elbowTh = 0
            self.shoulderV = 0
            self.elbowV = 0
            self.shoulderT = 0
            self.elbowT = 0
            self.handX = armlink * cos(self.shoulderTh) + forearmlink * cos(self.shoulderTh + self.elbowTh)
            self.handY = armlink * sin(self.shoulderTh) + forearmlink * sin(self.shoulderTh + self.elbowTh)


def genSuccessors(arm):
    """
    this will create a list of arm structs that can be acheived from the current state at this time step. I will create
    successors that will have at leat one of the acuators at a torque of 5 or -5  and the other with 11 other torques making up to
    29 possibilities in total.
    :param arm: this is an arm stuct of an arm that is in the inital position
    :return: a tuple of arm structs of possible successors after a timestep
    """
    succlist = []
    for eT in range(-5, 6, 1):
        for sT in(-5, 5):
            succlist.append(ArmStruct(arm, sT,eT))
    for sT in range(-4, 5, 1):
        for eT in (-5, 5):
            succlist.append(ArmStruct(arm, sT, eT))
    return tuple(filter(lambda x: x.shoulderV < 60*pi and x.elbowV < 60*pi, succlist))


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


def heu(handposition, goalstates):
    """

    :param handposition: this is a tuple of floats the (x,y, time) that describe the locaton of the hand in 2d*time space
    :param goalstates: this is a tuple of goal states of the same syntax as hand positions that are in the reach of the ball
    :return: an optomistic heuristic for the getting to the most optomistic goal state
    """
    best = float("inf")
    for goal in goalstates:
        if handposition[2] < goal[2]:
            best = min(best, (((handposition[0] - goal[0]) ** 2 + (handposition[1] - goal[1]) ** 2) ** 0.5) /(goal[2] - handposition[2]))

    return best


def bisectInsert(item, sortedList, key=lambda x: x):

    """
    this does insteriton of the item into the list to stay sorted
    :param item: item to insert into list
    :param sortedList: the list that will be inserted to
    :param key: this decides what to consider while sorting the list
    :return: a sorted list with the new item inserted
    """
    low = 0
    high = len(list)
    while True:
        if low == high:
            return sortedList[:low]+[item]+ sortedList[high:]
        ndx = int((high - low)/2)+ low
        if key(item) > key(sortedList[ndx]):
            low = ndx + 1
            continue
        high = ndx


def solutionSearch(initialArm, goalstates):
    """
    this will do a greedy search for a catch condition of the ball using a graph search
    :param initialArm: this i the initial state of the arm
    :param goalstates:
    :return:
    """
    fringe = [(0, initialArm, [])]
    while fringe:
        _, thisState, thispath = fringe.pop()
        if isGoal((thisState.handX, thisState.handY, thisState.time), goalstates):
            return thisState, thispath
        succs = genSuccessors(thisState)
        succs = [(heu(x), x, thispath+[x]) for x in succs]
        succs = list(filter(lambda x: x[0] != float("inf"), succs))
        #########THIS IS WHAT WILL TAKE MOST OF THE TIME############
        for item in succs:
            fringe = bisectInsert(item, fringe, key=lambda x: x[0])
        #########END################################################
    return None


def armmotion(vect, t):
    theta = vect[0]
    theta_dot = vect[1]
    tau = 5
    theta_ddot = (tau-(mass*g*radius*numpy.sin(theta))) / (i + pow(radius, 2))
    return [theta_dot, theta_ddot]


def main():
    while True:
        try:
            point1 = [float(x) for x in input("Enter first point in radial coordinates:" ).split(" ")]
            point2 = [float(x) for x in input("Enter second point in radial coordinates:").split(" ")]
            if len(point1) == 2 and len(point2) == 2:
                break
        except:print("Let's try this again you have to put in numeral numbers with or without a decimal")
    goals = findGoalsFromStart(point1, point2)
    if not goals:
        print("Ball. Does not reach the arms catch zone.")
        return None
    answer = catchsolutionSearch(ArmStruct(None, 0, 0), goals)
    if not answer:
        print("Strike. No catching motion was able to be identified.")
        return None
    catchstate, path = answer
    print("Out, for some reason hitting the ball is an out here like we are feilding even thouhg the rest of it is like batting.")
    print("the catch was done with arm postions ElbowTh:", catchstate.elbowTh, "ShoulderTh:", catchstate.shoulderTh,
          "and is a position that can be aceived as early as:", catchstate.time)


if __name__ == "__main__":
    main()