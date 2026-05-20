import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_csv("output/state_history.csv")

# 3D orbit plot
fig = plt.figure(figsize=(14, 4))

ax1 = fig.add_subplot(141, projection='3d')
ax1.plot(df.pos_x, df.pos_y, df.pos_z, lw=0.8)
ax1.set_title("Orbit (ECI)")
ax1.set_xlabel("X (m)"); ax1.set_ylabel("Y (m)"); ax1.set_zlabel("Z (m)")

# Velocity magnitude
ax2 = fig.add_subplot(142)
ax2.plot(df.vel_x, label="vx")
ax2.plot(df.vel_y, label="vy")
ax2.plot(df.vel_z, label="vz")
ax2.plot((df.vel_x**2 + df.vel_y**2 + df.vel_z**2)**0.5, label="|v|", lw=2, color='k')
ax2.legend()
ax2.set_title("Velocity")
ax2.set_xlabel("Timestep"); ax2.set_ylabel("m/s")

# Quaternion components
ax3 = fig.add_subplot(143)
ax3.plot(df.q_w, label="qw")
ax3.plot(df.q_x, label="qx")
ax3.plot(df.q_y, label="qy")
ax3.plot(df.q_z, label="qz")
ax3.set_title("Attitude (Quaternion)")
ax3.set_xlabel("Timestep"); ax3.set_ylabel("Unitless")
ax3.legend()

# Angular velocity
ax4 = fig.add_subplot(144)
ax4.plot(df.wx, label="ωx")
ax4.plot(df.wy, label="ωy")
ax4.plot(df.wz, label="ωz")
ax4.set_title("Angular Velocity")
ax4.set_xlabel("Timestep"); ax4.set_ylabel("rad/s")
ax4.legend()

plt.tight_layout()
plt.savefig("output/orbit_plot.png", dpi=150)
plt.show()

# Orbital Elements plot
df_oe = pd.read_csv("output/oe_history.csv")
angular_cols = df_oe.columns[2:]  # i, raan, argp, nu
df_oe[angular_cols] = df_oe[angular_cols] * (180 / np.pi)
fig, axs = plt.subplots(2, 3, figsize=(15, 8))
oe_labels = ["a (m)", "e", "i (deg)", "RAAN (deg)", "argp (deg)", "nu (deg)"]
for i, col in enumerate(df_oe.columns):
    ax = axs[i//3, i%3]
    ax.plot(df_oe[col])
    ax.set_title(oe_labels[i])
    ax.set_xlabel("Timestep")
plt.tight_layout()
plt.savefig("output/oe_plot.png", dpi=150)
plt.show()