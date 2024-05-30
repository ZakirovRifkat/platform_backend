from fillipov import fillipov
import numpy as np

from typing import List, Union
from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

origins = ["*"]  # Разрешить все источники

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/model", response_model=dict)
async def model(A: float, B: float, C: float, delta: float, initial: list = Query()):

    params = [A, B, C, delta]

    # Начальные условия
    y0 = list(map(float, initial))

    # Время интеграции
    T = 500
    t_span = [0, T]

    options = {"atol": 1e-10, "rtol": 1e-10}

    t, y, te, ye, ie, se = fillipov(
        vfields=vectorfields,
        jacobians=jacobians,
        pfunction=None,
        tspan=t_span,
        y0=y0,
        params=params,
        C=1,
        inopts=options,
    )
    t = t.tolist()
    y = y.tolist()
    result = {
        "t": t,
        "y": y,
    }

    return result


def jacobians(
    t,
    y,
    params,
):
    # Распаковка параметров
    A, B, C, delta = params

    # Якобиан в области S1 (H(x) > 0)
    J1 = np.array(
        [
            [0, 1, 0, 0],
            [-A, -B, 1, 0],
            [0, 0, -C, -1],
            [1 / delta, 0, 0, -1 / delta],
        ]
    )

    # Якобиан в области S2 (H(x) < 0)
    J2 = J1

    # Градиент градиента H, вектор, нормальный к поверхности разрыва
    d2H = np.zeros_like(J1)

    return J1, J2, d2H


def vectorfields(t, y, params):
    # Распаковка параметров
    A, B, C, delta = params

    # Векторное поле в области 1 - H(x) > 0
    F1 = np.array(
        [
            y[1],
            -B * y[1] - A * y[0] + y[2] - 1 / 2,
            -C * y[2] - y[3],
            (y[0] - y[3]) / delta,
        ]
    )

    # Векторное поле в области 2 - H(x) < 0
    F2 = np.array(
        [
            y[1],
            -B * y[1] - A * y[0] + y[2] + 1 / 2,
            -C * y[2] - y[3],
            (y[0] - y[3]) / delta,
        ]
    )
    # Переключающее многообразие
    H = y[1]

    # Вектор нормали к переключающему многообразию
    dH = np.array([0, 1, 0, 0])

    # Плоскость Пуанкаре
    h = 1

    # Направление расположения
    hdir = 1

    return F1, F2, H, dH, h, hdir
