import meep as mp
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd

# Parameters for the ring resonator
ring_outer_radius = 10  # outer radius of the ring
ring_width = 1  # width of the ring
ring_inner_radius = ring_outer_radius - ring_width  # inner radius of the ring
wg_width = 1  # width of the waveguide
gap = 0.1  # gap between waveguide and ring
dpml = 2  # PML thickness
cell_size = mp.Vector3(40, 40, 0)  # size of the simulation cell

# Medium definitions
ring_material = mp.Medium(epsilon=12)  # dielectric material for the ring
air = mp.Medium(epsilon=1)  # material for air (assumed to be surrounding medium)

# Geometry definitions
ring = mp.Cylinder(radius=ring_outer_radius, height=mp.inf, material=ring_material)
air_hole = mp.Cylinder(radius=ring_inner_radius, height=mp.inf, material=air)

# Waveguide definition
# wg = mp.Block(size=mp.Vector3(mp.inf, wg_width, mp.inf),
#               center=mp.Vector3(ring_inner_radius + gap + wg_width/2, 0),
#               material=ring_material)

wg = mp.Block(size=mp.Vector3(mp.inf, wg_width, mp.inf), 
              center=mp.Vector3(0, -ring_outer_radius - wg_width/2 - gap), 
              material=ring_material)

# PML layers
pml_layers = [mp.PML(dpml)]

# Source definition (continuous wave source for simplicity)
fcen = 0.1  # center frequency
df = 0.08

# Source definition (with an EigenModeSource for directionality)
source_position = -cell_size.x/2 + dpml + wg_width/2  # x position of the source
sources = [mp.EigenModeSource(src=mp.ContinuousSource(frequency=fcen),
                              size=mp.Vector3(0, wg_width, 0),
                              center=mp.Vector3(source_position, -ring_outer_radius - wg_width/2 - gap),
                              eig_band=1,
                              eig_parity=mp.ODD_Z+mp.EVEN_Y,
                              eig_match_freq=True)]

# Source definition (continuous wave source for simplicity)
# fcen = 0.1  # center frequency
# df = 0.08

# # Source definition (with an EigenModeSource for directionality)
# source_position = -cell_size.x/2 + dpml + wg_width/2  # x position of the source
# sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=df),
#                               size=mp.Vector3(0, wg_width, 0),
#                               center=mp.Vector3(source_position, -ring_outer_radius - wg_width/2 - gap),
#                               eig_band=1,
#                               eig_parity=mp.ODD_Z+mp.EVEN_Y,
#                               eig_match_freq=True)]

# Simulation
sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=[ring, wg, air_hole],
                    sources=sources,
                    resolution=10)

# Add flux regions
# input_flux_region = mp.FluxRegion(center=mp.Vector3(source_position, -ring_outer_radius - wg_width/2 - gap),
#                                   size=mp.Vector3(0, wg_width))
# output_flux_region = mp.FluxRegion(center=mp.Vector3(-source_position, -ring_outer_radius - wg_width/2 - gap),
#                                    size=mp.Vector3(0, wg_width))


# input_flux = sim.add_flux(fcen, df, 1, input_flux_region)
# output_flux = sim.add_flux(fcen, df, 1, output_flux_region)

# # Lists for time-dependent data
# times = []
# input_flux_values = []
# output_flux_values = []

# def record_flux(sim):
#     times.append(sim.meep_time())
#     input_flux_values.append(mp.get_fluxes(input_flux)[0])  # Assuming one frequency point
#     output_flux_values.append(mp.get_fluxes(output_flux)[0])

# sim.run(mp.at_every(1, record_flux), until=20000)  # Adjust the interval and duration as needed

# # Saving the flux data
# flux_df = pd.DataFrame({
#     'Time': times,
#     'Input Flux': input_flux_values,
#     'Output Flux': output_flux_values
# })
# flux_df.to_csv('flux_data.csv', index=False)


# Setup animation
animation = mp.Animate2D(sim, fields=mp.Ez, realtime=False, normalize=True)

# Run the simulation once and record all data
sim.run(mp.at_every(1, animation), until=2000)

# Generate and save the GIF
animation.to_gif(10, filename='simulation2.gif')

