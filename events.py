import numpy as np


def event0(t, y, vfields, jacobians, params, C, state, dir):

    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    value = H
    if state[0] == 1 or state[1] == 1:
        pass
    elif state[2] == 1:
        value = 1
    else:
        raise ("ERROR: Wrong event in filippov:event0")
    return value


def event1(t, y, vfields, jacobians, params, C, state, dir):
    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    dHF1 = dH @ F1

    value = dHF1
    if state[0] == 1 or state[1] == 1:
        pass
    elif state[2] == 1:
        pass
    else:
        raise ("ERROR: Wrong event in filippov:event1")

    return value


def event2(t, y, vfields, jacobians, params, C, state, dir):
    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    dHF2 = dH @ F2

    value = dHF2
    if state[0] == 1 or state[1] == 1:
        pass
    elif state[2] == 1:
        pass
    else:
        raise ("ERROR: Wrong event in filippov:event2")

    return value


def event3(t, y, vfields, jacobians, params, C, state, dir):
    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    value = h
    if state[0] == 1 or state[1] == 1:
        pass
    elif state[2] == 1:
        pass
    else:
        raise ("ERROR: Wrong event in filippov:event3")

    return value


def event4(t, y, vfields, jacobians, params, C, state, dir):
    F1, F2, H, dH, h, hdir = vfields(t, y, params)

    dHF1 = dH @ F1
    dHF2 = dH @ F2

    value = 1

    if state[0] == 1 or state[1] == 1:
        pass
    elif state[2] == 1:
        if jacobians is None:
            value = 1
        else:
            J1, J2, d2H = jacobians(t, y, params)

            dHF1_p_dHF2 = dHF1 + dHF2
            dHF2_dHF1 = dHF2 - dHF1

            Hu = -((dHF1_p_dHF2) / (dHF2_dHF1))

            F2_F1 = F2 - F1
            F2_F1_2 = F2_F1.dot(0.5)
            F1_p_F2 = F1 + F2
            F1_p_F2_2 = F1_p_F2.dot(0.5)

            J2_J1 = J2 - J1
            J2_J1_2 = J2_J1.dot(0.5)
            J1_p_J2 = J1 + J2
            J1_p_J2_2 = J1_p_J2.dot(0.5)

            # TODO Может косячу по формуле
            dHu = -(
                ((((F1_p_F2.T) @ d2H) + (dH @ J1_p_J2)) * dHF2_dHF1)
                - (((F2_F1.T) @ d2H) + (dH @ J2_J1)) * (dHF1_p_dHF2)
            ) / (dHF2_dHF1**2)

            F = (F1_p_F2_2) + (F2_F1_2) * Hu - dH.T.dot(C * H)

            dHuF = dHu @ F
            value = dHuF
    else:
        raise ("ERROR: Wrong event in filippov:event3")

    return value
