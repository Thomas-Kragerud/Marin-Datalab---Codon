



import time


# Defining variables
# ---------------------------------------------------------------
n = 30
epsi = 1e-6

nn = [0, 5, 10, 20, 30, 40, 60, 100, 500]
oo = [1.7, 1.78, 1.86, 1.92, 1.95, 1.96, 1.97, 1.98, 1.99]
#omega = scipy.interpolate.interp1d(nn, oo, kind='linear')(n)
omega = 1.95
Re = 100.0
tmax = 10.0
dt = 0.01
itmax = 30
h = 1 / n
beta = omega * (h ** 2) / (4 * dt)


u = [[0.0 for _ in range(n + 2)] for _ in range(n + 2)]
v = [[0.0 for _ in range(n + 2)] for _ in range(n + 2)]
p = [[0.0 for _ in range(n + 2)] for _ in range(n + 2)]




def frange(start, stop, step):
    current = start
    while current <= stop:
        yield current
        current += step

# ---------------------------------------------------------------
#@jit(nopython=True)

def funkesen(n, epsi, u : list[list[float]], v : list[list[float]], p, h, beta, Re, tmax, dt, itmax, omega, iter2=0, div=0, iflag=0) ->  list[list[float]]:

    for t in range(100):
        # Your loop body goes here

        #t = t + dt

        for i in range(1, n + 1):
            for j in range(1, n + 1):
                fux = ((u[i][j] + u[i + 1][j]) ** 2 - (u[i - 1][j] + u[i][j]) ** 2) * 0.25 / h

                fuy = ((v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1]) - (v[i][j - 1] + v[i + 1][j - 1]) * (
                            u[i][j - 1] + u[i][j])) * 0.25 / h

                fvx = ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i - 1][j] + u[i - 1][j + 1]) * (
                            v[i - 1][j] + v[i][j])) * 0.25 / h

                fvy = ((v[i][j] + v[i][j + 1]) ** 2 - (v[i][j - 1] + v[i][j]) ** 2) * 0.25 / h

                visu = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4.0 * u[i][j]) / (Re * h ** 2)
                visv = (v[i + 1][j] + v[i - 1][j] + v[i][j + 1] + v[i][j - 1] - 4.0 * v[i][j]) / (Re * h ** 2)

                u[i][j] = u[i][j] + dt * (((p[i][j] - p[i + 1][j]) / h) - fux - fuy + visu)
                v[i][j] = v[i][j] + dt * (((p[i][j] - p[i][j + 1]) / h) - fvx - fvy + visv)

        for it in range(itmax):
            iter2 = it

            for j in range(0, n + 2):
                u[0][j] = 0.0
                v[0][j] = - v[1][j]
                u[n][j] = 0.0
                v[n + 1][j] = - v[n][j]

            for i in range(0, n + 2):
                v[i][n] = 0.0
                v[i][0] = 0.0
                u[i][n + 1] = - u[i][n] + 2.0
                u[i][0] = - u[i][1]
            iflag = 0
            for it in range(0, itmax):
                iter2 = it

            for j in range(0, n + 2):
                u[0][j] = 0.0
                v[0][j] = - v[1][j]
                u[n][j] = 0.0
                v[n + 1][j] = - v[n][j]

            for i in range(0, n + 2):
                v[i][n] = 0.0
                v[i][0] = 0.0
                u[i][n + 1] = - u[i][n] + 2.0
                u[i][0] = - u[i][1]
            iflag = 0
            for j in range(1, n + 1):
                for i in range(1, n + 1):
                    div = (u[i][j] - u[i - 1][j]) / h + (v[i][j] - v[i][j - 1]) / h
                    if abs(div) >= epsi:
                        iflag = 1
                    delp = -beta * div
                    p[i][j] = p[i][j] + delp
                    u[i][j] = u[i][j] + delp * dt / h
                    u[i - 1][j] = u[i - 1][j] - delp * dt / h
                    v[i][j] = v[i][j] + delp * dt / h
                    v[i][j - 1] = v[i][j - 1] - delp * dt / h
            if iflag == 0:
                print("heo")
                break
        # if iter2 >= itmax:
        #     print('Warning', 'time: ', t, 'iter: ', iter2, 'Div: ', div)
        #
        # else:
        #     print('Time: ', t, 'iter: ', iter2)
    return u

start_time = time.time()
funkesen(n,epsi,u,v,p,h,beta,Re,tmax,dt,itmax,omega)
end_time = time.time()
print(end_time - start_time)