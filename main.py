from data import load_data
from utils import rl, c45, cart
import matplotlib.pyplot as plt


def methods(data):

    return results


data = load_data()

results = methods(data)

class_name = rl

if class_name == "rl":
    my_class = rl(results)
elif class_name == "c45":
    my_class = c45(results)
elif class_name == "cart":
    my_class = cart(results)
else:
    print("Invalid class name entered.")
    exit()

my_class.some_method()

results = my_class.run()

fig, axs = plt.subplots(nrows=2, ncols=2)

axs[0, 0].plot(results[0])
axs[0, 0].set_title("Plot 1")
axs[0, 1].plot(results[1])
axs[0, 1].set_title("Plot 2")
axs[1, 0].plot(results[2])
axs[1, 0].set_title("Plot 3")
axs[1, 1].plot(results[3])
axs[1, 1].set_title("Plot 4")

plt.show()