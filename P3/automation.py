import os

os.system("./twoBody.x 1 1000 50 1")
os.system("python3 plot.py twoBody Euler 50")

os.system("./twoBody.x 1 1000 100 1")
os.system("python3 plot.py twoBody Euler 100")

os.system("./twoBody.x 1 1000 200 1")
os.system("python3 plot.py twoBody Euler 200")

os.system("./twoBody.x 1 1000 400 1")
os.system("python3 plot.py twoBody Euler 400")


os.system("./twoBody.x 2 1000 50 1")
os.system("python3 plot.py twoBody Verlet 50")

os.system("./twoBody.x 2 1000 100 1")
os.system("python3 plot.py twoBody Verlet 100")

os.system("./twoBody.x 2 1000 200 1")
os.system("python3 plot.py twoBody Verlet 200")

os.system("./twoBody.x 2 1000 400 1")
os.system("python3 plot.py twoBody Verlet 400")
