import Tools

# start computation for *.off file, distribution ['cartesian', 'hexdens', 'cylindric'], particle radius
x = Tools.PlacePoints('/Vortex_Tube.off', 'cylindric', 4e-4)
# write ccordinates to csv
x.write_csv('points_inside.csv')