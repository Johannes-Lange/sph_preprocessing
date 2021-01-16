import Tools
import time
import pandas as pd
t1 = time.time()
x = Tools.PlacePoints('/Vortex_Tube.off', 'cartesian', 2e-4)
print(time.time()-t1)
"""
t1 = time.time()
x = Tools.PlacePoints('/Vortex_Tube.off', 'cartesian', 4e-4)
x.write_csv('points_inside_cart_4.csv')
t2 = time.time()

x = Tools.PlacePoints('/Vortex_Tube.off', 'cartesian', 2e-4)
x.write_csv('points_inside_cart_2.csv')
t3 = time.time()

x = Tools.PlacePoints('/Vortex_Tube.off', 'cartesian', 1e-4)
x.write_csv('points_inside_cart_1.csv')
t4 = time.time()

time = [t2-t1, t3-t2, t4-t3]
time = pd.DataFrame(time)
time.to_csv('cart.csv')
"""
