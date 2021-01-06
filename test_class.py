import Tools
import pandas as pd

x = Tools.ToBeNamed('/bunny.off', 'cylindric', 0.003)

out = [[p.x(), p.y(), p.z()] for p in x.points]
out = pd.DataFrame(out, columns=['x', 'y', 'z'])  # 45600

out.to_csv('points_inside_cylindric.csv')