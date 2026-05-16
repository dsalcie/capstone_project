def RK4(f, g, h, eos, t0, y0, z0) -> tuple:
    '''
    Implementation of Runge Kutta 4th Order-2nd ODE Solver, 1 iteration
    
    -Note: the TOV_solver function will take care of iterating through different r values.
     Runge Kutta 4th Order will take care of integrating for the given t
    
    -Takes:
        t: independent variable
        f: f(t, y, z) 
        g: g(t, y, z)
        h: step
        t0: initial value of independent variable
        y0: y(0), initial condition
        z0: z(0), initial condition
    -Returns: t, y, z values
    '''
    t_n = t0
    y_n = y0
    z_n = z0

    k_1y = f(t_n, y_n, z_n, eos)
    k_1z = g(t_n, y_n, z_n, eos)
    k_2y = f(t_n + h/2, y_n + ((h/2)*k_1y), z_n + ((h/2)*k_1z), eos)
    k_2z = g(t_n + h/2, y_n + ((h/2)*k_1y), z_n + ((h/2)*k_1z), eos)
    k_3y = f(t_n + h/2, y_n + ((h/2)*k_2y), z_n + ((h/2)*k_2z), eos)
    k_3z = g(t_n + h/2, y_n + ((h/2)*k_2y), z_n + ((h/2)*k_2z), eos)
    k_4y = f(t_n + h, y_n + (h * k_3y), z_n + (h * k_3z), eos)
    k_4z = g(t_n + h, y_n + (h * k_3y), z_n + (h * k_3z), eos)

    t_n += h
    y_n += (h/6)*(k_1y + (2*k_2y) + (2*k_3y) + k_4y)
    z_n += (h/6)*(k_1z + (2*k_2z) + (2*k_3z) + k_4z)
    return t_n, y_n, z_n