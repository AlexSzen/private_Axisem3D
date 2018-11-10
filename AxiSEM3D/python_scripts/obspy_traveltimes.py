from obspy.taup import TauPyModel

model = TauPyModel( model = "prem")
arrivals = model.get_travel_times(source_depth_in_km=55, distance_in_degree=70)

print(arrivals)