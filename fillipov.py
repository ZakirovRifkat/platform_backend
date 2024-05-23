import numpy as np
from scipy.integrate import solve_ivp
from events import event0,event1,event2,event3,event4


event0.terminal = 1
event1.terminal = 1
event2.terminal = 1
event3.terminal = 1
event4.terminal = 0

# event4.direction = 0

def get_events_number(yEvents):
    ie = np.array([], dtype=int)

    for i, arr in enumerate(yEvents):
        if(len(arr) !=0 ):
            ie = np.append(ie, i+1)
    return ie


def events_direction(t, y, vfields, jacobians, params, state, dir):
    F1, F2, H, dH, h, hdir = vfields(t, y, params)
    direction = np.copy(np.append(dir, hdir))
    event0.direction = direction[0]
    event1.direction = direction[1]
    event2.direction = direction[2]
    event3.direction = direction[3]
    
    if state[0] == 1 or state[1] == 1:
        direction[0] = -state[0]
        event0.direction = -state[0]
 
def fillipov(vfields, jacobians, pfunction, tspan, y0, params, C, inopts):
    t0 = tspan[0]
    t1 = tspan[-1]

    state, dir = findstate(vfields, jacobians, t0, y0, params)
    options = inopts

    yvect = np.array([])
    tvect = np.array([])
    te = np.array([])
    ye = np.array([])
    ie = np.array([])
    se = np.array([])
    
    stopit = False
    a = 0
    prev_IE = 0
    amount = 0
    while not stopit:
        
        # TODO: Вот тут вылезут ошибки, скорее всего
        
        sol = solve_ivp(filippovfunc, t_span=tspan, y0=y0, events=(event0, event1, event2, event3, event4), args=(vfields, jacobians, params, C, state, dir), method='DOP853', max_step=0.1, **options)

        # stopit = True
        a += 1
        
        t = sol.t
        y = sol.y.T
        TE = sol.t_events
        YE = sol.y_events
        IE = get_events_number(sol.t_events)
        if(len(IE) == 1 and prev_IE == IE[0]):
            amount +=1
        else:
            amount == 0
        if(len(IE)!=0):
            prev_IE = IE[0]
        if(amount > 100):
            raise("ERROR IN SYSTEM!!!")
        # print("t_events = ", TE)
        # print("y_events = ", YE)
        # print("IE = ", IE)
        # print("-----------------------------------\n")
        
        
        y0 = y[-1, :]

        # TODO: Здесь получается обычные списки питона
        yvect_list = list(yvect)
        for arr in y:
            yvect_list.append(arr)
        
        yvect = np.array(yvect_list)
        tvect = np.append(tvect, t)
        print("IE = ", IE)
        if len(IE) != 0:
            te = np.append(te, TE[IE[0]-1])
            ye = np.append(ye, YE[IE[0]-1])
            
        tspan = np.array([t[-1], t1])
        if len(IE) != 0 and (t[-1] != t1):
            for k in range(len(IE)):
                ie = np.append(ie, IE[k])
                if IE[k] == 4:
                    if pfunction:
                        y0 = pfunction(t, y0, params)
                else:
                    if state[2] == 1:
                        
                        if IE[k] == 2 or IE[k] == 3:
                            state[IE[k] - 2] = -state[IE[k] - 2]
                            state[2] = -state[2]
                            state[3] = -state[3]
                            state[4] = -state[4]
                            
                            dir[[0, IE[k] - 1]] =  np.array([-1, -dir[IE[k] - 1]])
                            
                        elif IE[k] == 5:
                            pass  
                        else:
                            raise('Error, there is something wrong with the event in filippov')
                    elif state[3] == 1:
                        
                        if IE[k] == 1:
                            state[0] = -state[0]
                            state[1] = -state[1]
                            dir[IE[k] - 1] = -dir[IE[k] - 1]
                        elif IE[k] == 2 or IE[k] == 3:
                            state[3] = -state[3]
                            state[4] = -state[4]
                            dir[IE[k] - 1] = -dir[IE[k] - 1]
                        elif IE[k] == 5:
                            pass 
                        else:
                            raise('Error, there is something wrong with the event in filippov')
                    elif state[4] == 1:
                        
                        if IE[k] == 1:
                            state[0] = -1
                            state[1] = -1
                            state[2] = -state[2]
                            dir[IE[k] - 1 ] = -dir[IE[k] - 1]
                            
                        elif IE[k] == 2 or IE[k] == 3:
                            state[3] = -state[3]
                            state[4] = -state[4]
                            dir[IE[k] - 1] = -dir[IE[k] - 1]
                        elif IE[k] == 5:
                            pass 
                        else:
                            raise("Error, there is something wrong with the event in filippov")
                    else:
                        raise("Error, There is something wrong with the state vector in filippov")
                se = np.append(se, state)
        elif IE:
            stopit = True
            ie = np.append(ie, IE)
            se = np.append(se, state)
        else:
            stopit = True
                        
                
    return tvect, yvect, te, ye, ie, se



# ---------------------------------------------------------------

def findstate(vfields, jacobians, t0, y0, params):
    state = np.full(5, -1)

    """
        vfields: Векторное поле
        F1 - вектор
        F2 - вектор 
        H - переменная по которой идет разрыв
        dH - вектор нормали к переменной разрыва
        h1 - число = 1
        hdir - число = 1 
    """
    F1, F2, H, dH, h1, hdir = vfields(t0, y0, params)

    dHF1 = np.dot(dH, F1)
    dHF2 = np.dot(dH, F2)

    dir = np.array([-np.sign(H), -np.sign(dHF1), -np.sign(dHF2)])

    if H > 0:
        state[0] = -state[0]
    elif H < 0:
        state[1] = -state[1]
    elif np.sign(dHF1) * np.sign(dHF2) < 0:
        state[2] = -state[2]
    else:
        if np.sign(dHF1) > 0:
            state[0] = -state[0]
        else:
            state[1] = -state[1]

    if np.sign(dHF1) * np.sign(dHF2) > 0:
        state[3] = -state[3]
    elif np.sign(dHF1) * np.sign(dHF2) < 0:
        state[4] = -state[4]
    else:
        if jacobians is None:
            state[3] = -state[3]
        else:
            """
                Тут все вектора
            """
            J1, J2, d2H = jacobians(t0, y0, params)

            if dHF1 == 0:
                HxF1x_F1Hxx = dH @ J1 + F1.T @ d2H
                sig = np.sign(HxF1x_F1Hxx @ F1) * np.sign(dHF2)
                dir[1] = -np.sign(HxF1x_F1Hxx @ F1)
            elif dHF2 == 0:
                HxF2x_F2Hxx = dH @ J2 + F2.T @ d2H
                sig = np.sign(HxF2x_F2Hxx @ F2) * np.sign(dHF1)
                dir[2] = -np.sign(HxF2x_F2Hxx @ F2)
            else:
                sig = 1
                raise('ERROR: Something is wrong in filippov:findstate')
            
            if sig < 0:
                state[4]= -state[4]
            else:
                state[3]= -state[3]
    return state, dir



# --------------------------------------------------------------

def fevents(t, y, vfields, jacobians, params, C, state, dir):
    """
        vfields: Векторное поле
        F1 - вектор
        F2 - вектор 
        H - переменная по которой идет разрыв
        dH - вектор нормали к переменной разрыва
        h1 - число = 1
        hdir - число = 1 
    """
    print("t = ", t)
    print("y = ", y)

    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    dHF1 = dH @ F1
    dHF2 = dH @ F2
    

    value = np.array([H, dHF1, dHF2, h])
    direction = np.copy(np.append(dir, hdir))
    
    if state[0] == 1 or state[1] == 1: 
        direction[0] = -state[0]
    elif state[2] == 1:
        value[0] = 1
        if jacobians is None:
            value = np.append(value, 1)
            direction = np.append(direction, 0)
        else:
            """
                Тут все вектора
            """
            J1, J2, d2H = jacobians(t, y, params)
            
            dHF1_p_dHF2 = dHF1 + dHF2
            dHF2_dHF1 = dHF2 - dHF1
            
            Hu =- ((dHF1_p_dHF2) / (dHF2_dHF1))
            
            F2_F1 = F2 - F1
            F2_F1_2 = 0.5 * F2_F1
            F1_p_F2 = F1 + F2
            F1_p_F2_2 = 0.5 * F1_p_F2

            J2_J1 = J2 - J1
            J2_J1_2 = 0.5 * J2_J1
            J1_p_J2 = J1 + J2
            J1_p_J2_2 = 0.5 * J1_p_J2
            
            
            # TODO Может косячу по формуле
            dHu = -((((F1_p_F2.T) @ d2H) + (dH @ J1_p_J2) @ dHF2_dHF1) - \
                    ((F2_F1.T @ d2H) + (dH @ J2_J1)) @ (dHF1_p_dHF2)) / (dHF2_dHF1 ** 2)

            F = (F1_p_F2_2) + (F2_F1_2) @ Hu - C @ H @ dH.T
            
            dHuF = (dHu @ F)
            value = np.append(value, dHuF)
            direction = np.append(direction, 0)
    else:
        raise("ERROR: Wrong event in filippov:fevents")
            
    isterminal = [1, 1, 1, 1, 0]
    
    print("------------")
    print("value = ", value)
    print("------------")
    return value, isterminal, direction


# -------------------------------------------------------------------

def filippovfunc(t, y, vfields, jacobians, params, C, state, dir):
    dy = np.zeros(len(y))
    """
        vfields: Векторное поле
        F1 - вектор
        F2 - вектор 
        H - переменная по которой идет разрыв
        dH - вектор нормали к переменной разрыва
        h1 - число = 1
        hdir - число = 1 
    """
    F1, F2, H, dH, h1, hdir = vfields(t, y, params)
    
    events_direction(t,y,vfields, jacobians, params, state, dir)
    if state[0] == 1:
        F = F1
    elif state[1] == 1:
        F = F2
    elif state[2] == 1:
        Fa = F1.dot(0.5)
        Fb = F2.dot(0.5)
        dHF1 = np.dot(dH,F1)
        dHF2 = np.dot(dH, F2)
        Hu =- ((dHF1 + dHF2) / (dHF2 - dHF1))
        F = (Fa + Fb) + Hu * (Fb - Fa) - C * H * dH.T
    else:
        raise("Error, there is something wrong with the state vector in filippov:filippovfunc") 
    dy = F 
    return dy
