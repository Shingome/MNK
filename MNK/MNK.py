import pandas as pd
import sympy as sp
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend


def get_sympy_subplots(plot):
    backend = MatplotlibBackend(plot)
    backend.process_series()
    backend.fig.tight_layout()
    return backend.ax[0]


def from_fun(fun, table):
    return list(fun.subs("x", i) for i in table["x"])


def calculate_e(y, fun, size):
    return [(y[i] - fun[i]) ** 2 for i in range(size)] + [None]


def plot(fun, x_list, y_list):
    pl = sp.plotting.plot(fun, show=False)
    axe = get_sympy_subplots(pl)
    axe.plot(x_list, y_list, 'o')
    return pl


def to_float(list):
    return (sp.Float(i) for i in list)


# Изначальная таблица
table = pd.DataFrame({"x": list(range(0, 11, 2)),
                      "y": [5, -1, 0.5, 1.5, 4.5, 8.5]})
n = len(table["x"])
x_list = (table['x'])
y_list = (table['y'])
print(table)

# Точечный график
plt.scatter(table["x"], table["y"])

# Дополнение таблицы для поиска параметров
for i in range(2, 7):
    table["x^{0}".format(i)] = table['x'] ** i
table["x*y"] = table["x"] * table["y"]
table["x^2*y"] = table["x^2"] * table["y"]
table["x^3*y"] = table["x^3"] * table["y"]
table.loc[len(table)] = list(table[i].sum() for i in table.columns)
table = table.rename(index={len(table) - 1: "sum"})
print(table)

# Инициализация переменных сигм для понятного представления математической формулы
sum_x2 = table["x^2"]["sum"]
sum_x3 = table["x^3"]["sum"]
sum_x4 = table["x^4"]["sum"]
sum_x5 = table["x^5"]["sum"]
sum_x6 = table["x^6"]["sum"]
sum_xy = table["x*y"]["sum"]
sum_x2y = table["x^2*y"]["sum"]
sum_x3y = table["x^3*y"]["sum"]
sum_x = table["x"]["sum"]
sum_y = table["y"]["sum"]

# Решение системы уравнений
a, b = sp.symbols('a, b')
solve = sp.linsolve([a * sum_x2 + b * sum_x - sum_xy,
                     a * sum_x + b * n - sum_y],
                    (a, b))
a, b = to_float(solve.args[0])

# Первый график
x = sp.symbols('x')
fun_1 = a * x + b
pl_1 = plot(fun_1, x_list, y_list)

# Решение второй системы уравнений
a, b, c = sp.symbols('a, b, c')
solve = sp.linsolve([a * sum_x4 + b * sum_x3 + c * sum_x2 - sum_x2y,
                     a * sum_x3 + b * sum_x2 + c * sum_x - sum_xy,
                     a * sum_x2 + b * sum_x + c * n - sum_y],
                    (a, b, c))
a, b, c = to_float(solve.args[0])

# Второй график
fun_2 = a * x * x + b * x + c
pl_2 = plot(fun_2, x_list, y_list)

# Решение третьей системы уравнений
a, b, c, d = sp.symbols('a, b, c, d')
solve = sp.linsolve([a * sum_x6 + b * sum_x5 + c * sum_x4 + d * sum_x3 - sum_x3y,
                     a * sum_x5 + b * sum_x4 + c * sum_x3 + d * sum_x2 - sum_x2y,
                     a * sum_x4 + b * sum_x3 + c * sum_x2 + d * sum_x - sum_xy,
                     a * sum_x3 + b * sum_x2 + c * sum_x + d * n - sum_y],
                    (a, b, c, d))
a, b, c, d = to_float(solve.args[0])

# Третий график
fun_3 = (a * x ** 3) + (b * x ** 2) + (c * x) + d
pl_3 = plot(fun_3, x_list, y_list)

# Вычисление квадратов отклонений
table = table[["x", "y"]][:len(table['x']) - 1]
table.loc[len(table)] = list(None for i in table.columns)
table = table.rename(index={len(table) - 1: "Q"})
for i in ((1, fun_1), (2, fun_2), (3, fun_3)):
    table["f{0}(x)".format(i[0])] = from_fun(i[1], table)
for i in range(1, 4):
    e, f = "e{0}".format(i), 'f{0}(x)'.format(i)
    table[e] = calculate_e(table['y'], table[f], n)
    table.at["Q", e] = table[e].sum()

# Совмещение графиков
final_plot = pl_1
final_plot.append(pl_2[0])
final_plot.append(pl_3[0])
axe = get_sympy_subplots(final_plot)
axe.plot(x_list, y_list, 'o')
plt.show()
