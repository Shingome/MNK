import pandas as pd
import sympy as sp
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot


def get_sympy_subplots(plot: Plot):
    backend = MatplotlibBackend(plot)

    backend.process_series()
    backend.fig.tight_layout()
    return backend.ax[0]


# Изначальная таблица
table = pd.DataFrame({"x": list(range(0, 11, 2)),
                      "y": [5, -1, 0.5, 1.5, 4.5, 8.5]})
n = len(table["x"])
print(table)

# Точечный график
plt.scatter(table["x"], table["y"])

# Дополнение таблицы для поиска параметров
table["x^2"] = table["x"] ** 2
table["x^3"] = table["x"] ** 3
table["x^4"] = table["x"] ** 4
table["x^5"] = table["x"] ** 5
table["x^6"] = table["x"] ** 6
table["x*y"] = table["x"] * table["y"]
table["x^2*y"] = table["x^2"] * table["y"]
table["x^3*y"] = table["x^3"] * table["y"]
table.loc[len(table)] = list(table[i].sum() for i in table.columns)
table = table.rename(index={len(table) - 1: "sum"})
print(table)

# Инициализация переменных сигм для облегчения записи системы
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
a, b = solve.args[0]
a, b = sp.Float(a), sp.Float(b)

# Первый график
x = sp.symbols('x')
fun_1 = a * x + b
pl_1 = sp.plotting.plot(fun_1, show=False)
axe = get_sympy_subplots(pl_1)
x_list = list(table["x"][:len(table['x']) - 1])
y_list = list(table["y"][:len(table['x']) - 1])
axe.plot(x_list, y_list, 'o')

# Решение второй системы уравнений
a, b, c = sp.symbols('a, b, c')
solve = sp.linsolve([a * sum_x4 + b * sum_x3 + c * sum_x2 - sum_x2y,
                     a * sum_x3 + b * sum_x2 + c * sum_x - sum_xy,
                     a * sum_x2 + b * sum_x + c * n - sum_y],
                    (a, b, c))
a, b, c = solve.args[0]
a, b, c = sp.Float(a), sp.Float(b), sp.Float(c)

# Второй график
fun_2 = a * x * x + b * x + c
pl_2 = sp.plotting.plot(fun_2, show=False)
axe = get_sympy_subplots(pl_2)
axe.plot(x_list, y_list, 'o')

# Решение третьей системы уравнений
a, b, c, d = sp.symbols('a, b, c, d')
solve = sp.linsolve([a * sum_x6 + b * sum_x5 + c * sum_x4 + d * sum_x3 - sum_x3y,
                     a * sum_x5 + b * sum_x4 + c * sum_x3 + d * sum_x2 - sum_x2y,
                     a * sum_x4 + b * sum_x3 + c * sum_x2 + d * sum_x - sum_xy,
                     a * sum_x3 + b * sum_x2 + c * sum_x + d * n - sum_y],
                    (a, b, c, d))
a, b, c, d = solve.args[0]
a, b, c, d = sp.Float(a), sp.Float(b), sp.Float(c), sp.Float(d)

# Третий график
fun_3 = (a * x ** 3) + (b * x ** 2) + (c * x) + d
pl_3 = sp.plotting.plot(fun_3, show=False)
axe = get_sympy_subplots(pl_3)
axe.plot(x_list, y_list, 'o')

# Вычисление квадратов отклонений
table_2 = table[["x", "y"]][:len(table['x']) - 1]
table_2.loc[len(table_2)] = list(None for i in table_2.columns)
table_2 = table_2.rename(index={len(table_2) - 1: "Q"})
table_2["f1(x)"] = list(fun_1.subs("x", i) for i in table_2["x"])
table_2["f2(x)"] = list(fun_2.subs("x", i) for i in table_2["x"])
table_2["f3(x)"] = list(fun_3.subs("x", i) for i in table_2["x"])
table_2["e1"] = list((table_2["y"][i] - table_2["f1(x)"][i]) ** 2 for i in range(n)) + [None]
table_2["e1"]["Q"] = table_2["e1"][:len(table['x']) - 1].sum()
table_2["e2"] = list((table_2["y"][i] - table_2["f2(x)"][i]) ** 2 for i in range(n)) + [None]
table_2["e2"]["Q"] = table_2["e2"][:len(table['x']) - 1].sum()
table_2["e3"] = list((table_2["y"][i] - table_2["f3(x)"][i]) ** 2 for i in range(n)) + [None]
table_2["e3"]["Q"] = table_2["e3"][:len(table['x']) - 1].sum()
print(table_2)

exit()

# Совмещение графиков
final_plot = pl_1
final_plot.append(pl_2[0])
final_plot.append(pl_3[0])
axe = get_sympy_subplots(final_plot)
axe.plot(x_list, y_list, 'o')
plt.show()