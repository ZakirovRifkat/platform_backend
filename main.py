import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from scipy.integrate import solve_ivp
from fillipov import fillipov


app = FastAPI()

origins = ["*"]  # Разрешить все источники

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/servomodel", response_model=dict)
async def servomodel(
    A: float,
    B: float,
    C: float,
    delta: float,
    initial: list[float] = Query(..., min_length=4, max_length=4),
):
    if delta == 0:
        raise HTTPException(status_code=422, detail="delta must not be 0")

    scalar_values = {
        "A": A,
        "B": B,
        "C": C,
        "delta": delta,
    }
    for name, value in scalar_values.items():
        if not np.isfinite(value):
            raise HTTPException(status_code=422, detail=f"{name} must be finite")

    if not all(np.isfinite(value) for value in initial):
        raise HTTPException(
            status_code=422, detail="initial must contain finite values"
        )

    params = [A, B, C, delta]

    # Начальные условия
    y0 = [float(value) for value in initial]

    # Время интеграции
    T = 500
    t_span = [0, T]

    options = {"atol": 1e-10, "rtol": 1e-10}

    try:
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
    except Exception as exc:
        raise HTTPException(
            status_code=500, detail=f"model calculation failed: {exc}"
        ) from exc

    t = t.tolist()
    y = y.tolist()
    result = {
        "t": t,
        "y": y,
    }

    return result


def watt_governor(variables, a, b, F0, m, J, beta, r, gamma0, x0):
    omega, y, z = variables

    return np.array(
        [
            y,
            z,
            -a * z - b * y - (F0 / (m * J)) * (beta * m * r * omega**2 - gamma0 * x0),
        ],
        dtype=float,
    )


@app.get("/classicmodel", response_model=dict)
async def classicmodel(
    a: float,
    b: float,
    F0: float,
    m: float,
    J: float,
    beta: float,
    r: float,
    gamma0: float,
    x0: float,
    dt: float = Query(0.01, gt=0),
    tEnd: float = Query(99.99, gt=0),
    initial: list[float] = Query(..., min_length=3, max_length=3),
):
    if m == 0:
        raise HTTPException(status_code=422, detail="m must not be 0")
    if J == 0:
        raise HTTPException(status_code=422, detail="J must not be 0")

    scalar_values = {
        "a": a,
        "b": b,
        "F0": F0,
        "m": m,
        "J": J,
        "beta": beta,
        "r": r,
        "gamma0": gamma0,
        "x0": x0,
        "dt": dt,
        "tEnd": tEnd,
    }
    for name, value in scalar_values.items():
        if not np.isfinite(value):
            raise HTTPException(status_code=422, detail=f"{name} must be finite")

    if not all(np.isfinite(value) for value in initial):
        raise HTTPException(
            status_code=422, detail="initial must contain finite values"
        )

    t_eval = np.arange(0.0, tEnd + dt * 0.5, dt, dtype=float)
    y0 = np.array([float(value) for value in initial], dtype=float)

    def rhs(_t, y):
        return watt_governor(y, a, b, F0, m, J, beta, r, gamma0, x0)

    try:
        solution = solve_ivp(
            fun=rhs,
            t_span=(0.0, float(tEnd)),
            y0=y0,
            t_eval=t_eval,
            method="RK45",
            rtol=1e-9,
            atol=1e-9,
            vectorized=False,
        )
    except Exception as exc:
        raise HTTPException(
            status_code=500, detail=f"classic model calculation failed: {exc}"
        ) from exc

    if not solution.success:
        raise HTTPException(
            status_code=500,
            detail=f"classic model solver failed: {solution.message}",
        )

    y_scipy = solution.y.T

    return {
        "t": solution.t.tolist(),
        "y": y_scipy.tolist(),
    }


def jacobians(
    _t,
    _y,
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


@app.get("/health", response_model=dict)
async def health():
    return {"status": "ok"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=False)
