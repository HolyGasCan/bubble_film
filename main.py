from matplotlib import pyplot as plt
import numpy as np
import time
from screeninfo import get_monitors
import copy

in_to_mm = 25.4
window_coeff = 0.8

screen = [m for m in get_monitors() if m.is_primary][0]
window_width = screen.width_mm / in_to_mm * window_coeff
window_height = screen.height_mm / in_to_mm * window_coeff

n = 100                 # amount of dots
height = 4              # total height of cylinder
r = 4                   # initial radius
dh = height / n         # delta height
dr = 0.001              # delta radius
eps = 1 / (10 ** 10)    # error
bubble = [[],           # r
          []]           # height


# theoretical curve
min_theor = 3.39335
x = np.linspace(0, 4)
y = min_theor * np.cosh((x - height / 2) / min_theor)


# initiating cylinder
for i in range(n + 1):
    bubble[0].append(r)
    bubble[1].append(i * dh)


# surface of a cylinder segment
def cyl_surface(r1, r2):
    return (r1 + r2) * np.pi * np.sqrt(1 + ((r2 - r1) / dh) ** 2) * dh


# calc overall surface with changed radius
def surface(func, index, r1):
    sur = 0
    cyl = copy.deepcopy(func)
    cyl[0][index] = r1
    cyl[0] = cyl[0][0:index] + [r1] * (len(cyl[0]) - 2 * index) + cyl[0][len(cyl[0])-index:len(cyl[0])]

    for i in range(len(cyl[0]) - 1):
        sur += cyl_surface(cyl[0][i], cyl[0][i + 1])

    return sur


# change r to minimize surface
def minimize(func, index):
    r1 = func[0][index]
    m_dr = surface(func, index, r1 - dr)
    p_dr = surface(func, index, r1 + dr)

    if abs(m_dr - p_dr) < eps or r1 - dr < 0:
        return "min"

    current_surface = surface(func, index, r1)

    if m_dr < current_surface:
        return r1 - dr
    elif p_dr < current_surface:
        return r1 + dr
    else:
        return "min"


if __name__ == "__main__":
    minimized = ["no"] * (int(n / 2))
    isPopped = False
    cnt = 1

    plt.ion()
    fig = plt.figure(figsize=(window_width, window_height))
    graph = fig.add_subplot(111)
    max_x = max(bubble[0][0], bubble[0][n])
    graph.axis([-max_x - 1, max_x + 1, -1, height + 1])
    graph.axvline(ymin=1/(height+2), ymax=1-1/(height+2), ls='--', color='cyan')

    line_bot, = graph.plot([-bubble[0][0], bubble[0][0]], [0, 0],
                           color='blue')
    line_top, = graph.plot([-bubble[0][n], bubble[0][n]], [height, height],
                           color='blue')
    curve_right, = graph.plot(bubble[0], bubble[1],
                              color='blue')
    curve_left, = graph.plot([-x for x in bubble[0]], bubble[1],
                             color='blue')
    cosh_right, = graph.plot(y, x,
                             c='red')

    while "no" in minimized:
        for i in range(1, int(n / 2) + 1):
            resp = minimize(bubble, i)
            if resp == "min":
                minimized[i - 1] = "min"
            else:
                bubble[0][i] = resp
                bubble[0][n - i] = resp
                minimized[i - 1] = "no"
        if min(bubble[0]) < 2 * eps:
            isPopped = True

        curve_right.set_xdata(bubble[0])
        curve_left.set_xdata([-x for x in bubble[0]])

        if isPopped:
            graph.lines.pop(3)
            graph.lines.pop(3)
            fig.canvas.draw()
            fig.canvas.flush_events()
            graph.set_title("Bubble popped!")
            break

        graph.set_title("Calculating" + "." * (cnt // 5 % 4))
        fig.canvas.draw()
        fig.canvas.flush_events()
        cnt += 1
        time.sleep(0.01)

    if not isPopped:
        graph.set_title("Finished!")

    print(abs(min(bubble[0]) - min_theor))
    plt.ioff()
    plt.show()
