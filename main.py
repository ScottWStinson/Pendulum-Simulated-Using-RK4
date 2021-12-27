# This should simulate the motion of a simple pendulum
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sin, cos, pi


class Pendulum:

    def __init__(self, mass, gravity, length, time, theta_0, omega0, step):
        self.mass = mass
        self.gravity = gravity
        self.length = length
        self.g_over_L = -gravity / length
        self.h = step
        self.time = time
        self.theta0 = theta_0
        self.posx1 = []
        self.posy1 = []
        self.posx2 = []
        self.posy2 = []
        self.omega = [omega0]
        self.theta1 = [theta_0]
        self.theta2 = [theta_0]

    def start(self):
        self.non_linear_pendulum()
        self.linear_pendulum()
        self.linear_animator()
        self.non_linear_animator()

        return self.theta1, self.theta2

    def stop(self):
        self.ani1.event_source.stop()
        self.ani2.event_source.stop()

    def linear_pendulum(self):

        for i, v in enumerate(self.time):
            self.theta1.append(self.theta0 * cos(self.time[i]/10 * (self.gravity / self.length) ** 0.5))
            self.posx1.append(sin(self.theta1[i]) * self.length)
            self.posy1.append(cos(self.theta1[i]) * -self.length)

        return self.posx1, self.posy1, self.theta1, self.time

    def non_linear_pendulum(self):
        for k, v in enumerate(self.time):
            m1 = self.omega[k]
            k1 = self.g_over_L * sin(self.theta2[k])

            m2 = self.omega[k] + (0.5 * self.h * k1)
            k2 = self.g_over_L * sin(self.theta2[k] + (0.5 * self.h * m1))

            m3 = self.omega[k] + (0.5 * self.h * k2)
            k3 = self.g_over_L * sin(self.theta2[k] + (0.5 * self.h * m2))

            m4 = self.omega[k] + (self.h * k3)
            k4 = self.g_over_L * sin(self.theta2[k] + (self.h * m3))

            self.theta2.append(self.theta2[k] + ((self.h / 6) * (m1 + 2 * m2 + 2 * m3 + m4)))
            self.omega.append(self.omega[k] + ((self.h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)))
            self.posx2.append(sin(self.theta2[k]) * L)
            self.posy2.append(cos(self.theta2[k]) * -L)
            print(self.posx2[k], self.posy2[k])

        return self.posx2, self.posy2, self.theta2, self.time

    def plotter1(self, k):

        tail_length = 10
        plt.subplot(211)
        plt.cla()
        plt.title('Pendulum Simulation Using Small Angle Approximation')

        if k <= tail_length:
            for j in range(tail_length):
                if k <= j:
                    plt.plot(self.posx1[k], self.posy1[k], 'o', color='b', alpha=1 - (k / tail_length))
                    plt.plot([0, self.posx1[k]], [0, self.posy1[k]],
                             color='r', marker='o', linestyle='-')
                else:
                    plt.plot(self.posx1[k - j], self.posy1[k - j], 'o', color='b', alpha=1 - (j / tail_length))
                    plt.plot([0, self.posx1[k]], [0, self.posy1[k]],
                             color='r', marker='o', linestyle='-')

        if k > tail_length:
            for j in range(tail_length):
                plt.plot(self.posx1[k - j], self.posy1[k - j], 'o', color='b', alpha=1 - (j / tail_length))
                plt.plot([0, self.posx1[k]], [0, self.posy1[k]],
                         color='r', marker='o', linestyle='-')

        plt.xlim([-self.length * 1.5, self.length * 1.5])
        plt.ylim([-self.length * 1.5, self.length])
        plt.ylabel("Meters(m)")
        plt.grid(True)
        plt.tight_layout()

        if k >= len(self.time) - 1:
            self.stop()

    def plotter2(self, k):

        tail_length = 10
        plt.subplot(212)
        plt.cla()
        plt.title('Pendulum Simulation Using RK4')

        if k <= tail_length:
            for j in range(tail_length):
                if k <= j:
                    plt.plot(self.posx2[k], self.posy2[k], 'o', color='b', alpha=1 - (k / tail_length))
                    plt.plot([0, self.posx2[k]], [0, self.posy2[k]],
                             color='r', marker='o', linestyle='-')
                else:
                    plt.plot(self.posx2[k - j], self.posy2[k - j], 'o', color='b', alpha=1 - (j / tail_length))
                    plt.plot([0, self.posx2[k]], [0, self.posy2[k]],
                             color='r', marker='o', linestyle='-')
        if k > tail_length:
            for j in range(tail_length):
                plt.plot(self.posx2[k - j], self.posy2[k - j], 'o', color='b', alpha=1 - (j / tail_length))
                plt.plot([0, self.posx2[k]], [0, self.posy2[k]],
                         color='r', marker='o', linestyle='-')

        plt.xlim([-self.length * 1.5, self.length * 1.5])
        plt.ylim([-self.length * 1.5, self.length])
        plt.xlabel("Meters (m)")
        plt.ylabel("Meters (m)")
        plt.grid(True)
        plt.tight_layout()

        if k >= len(self.time) - 1:
            self.stop()

    def linear_animator(self):
        self.ani1 = FuncAnimation(plt.gcf(), self.plotter1, interval=100)

    def non_linear_animator(self):
        self.ani2 = FuncAnimation(plt.gcf(), self.plotter2, interval=100)


m = 1   # mass = 1 kg
g = 9.81    # gravity = 9.81 m/s^2
L = 1   # Length of string = 2 m
totalT = 50
t = [i for i in range(totalT)]
theta0 = pi/4
h = 0.1
omega0 = 0

attempt = Pendulum(m, g, L, t, theta0, omega0, h)
plt.figure(1)
SAT, RK4T = attempt.start()     #SAT = Small Angle Theta, RK4T is 4th order Runge Kutta Theta
plt.show()
# Comparing the angles for each case
plt.figure(2)
SAT.pop()
RK4T.pop()
plt.plot(t, SAT, color='m', linestyle='dotted')
plt.plot(t, RK4T, color='g', linestyle='dotted')
plt.xlabel('Time (deca-seconds)')
plt.ylabel('Angle from Vertical (rad)')
plt.grid(True)
plt.legend(['Small Angle Approximation', '4th Order Runge Kutta Approximation'], loc='lower right')
plt.show()

