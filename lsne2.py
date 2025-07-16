import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import time

# --- Constants (can be adjusted for visual effect) ---
# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11
# Speed of light (m/s)
c = 2.998e8

# Conversion factors for astronomical units to SI units
M_SUN = 1.989e30  # Mass of the Sun in kg
R_SUN = 6.957e8   # Radius of the Sun in meters
M_JUPITER = 1.898e27 # Mass of Jupiter in kg
R_JUPITER = 6.9911e7 # Radius of Jupiter in meters

# --- Page Configuration and Styling ---
st.set_page_config(
    page_title="Gravitational Lensing Simulator",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for universe theme and element styling
st.markdown(
    """
    <style>
    /* Main app background with a plain black color */
    .stApp {
        background-color: black; /* Changed to plain black */
        color: white; /* Default text color */
    }
    /* Sidebar background with some transparency */
    .css-1d391kg {
        background-color: rgba(0, 0, 0, 0.7);
    }
    /* Main content area background with some transparency and rounded corners */
    .css-1lcbmhc {
        background-color: rgba(0, 0, 0, 0.5);
        padding: 20px;
        border-radius: 10px;
    }
    /* Ensure all text elements are white */
    h1, h2, h3, h4, h5, h6, .stMarkdown, label, .stSlider > div > div > div {
        color: white !important;
    }
    /* Style for the slider fill */
    .stSlider > div > div > div > div {
        background-color: #4CAF50; /* Green color for the slider track */
    }
    /* Style for the start simulation button */
    .stButton > button {
        background-color: #4CAF50; /* Green background */
        color: white; /* White text */
        border-radius: 10px; /* Rounded corners */
        padding: 10px 20px; /* Padding inside the button */
        font-size: 1.2em; /* Larger font size */
        border: none; /* No border */
        cursor: pointer; /* Pointer on hover */
        transition: background-color 0.3s ease; /* Smooth transition on hover */
    }
    .stButton > button:hover {
        background-color: #45a049; /* Darker green on hover */
    }
    /* Style for the checkbox label */
    .stCheckbox > label {
        color: white;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🌌 Gravitational Lensing Simulator")
st.markdown("Explore how massive objects bend light and magnify distant stars.")

st.sidebar.header("Simulation Parameters")

# --- Sliders and Checkbox for user input ---

# Slider for the mass of the lensing star (Star 2)
mass_star2_solar = st.sidebar.slider(
    "Mass of Lensing Star (Star 2) (Solar Masses)",
    min_value=0.1, max_value=10.0, value=1.0, step=0.1,
    help="Adjust the mass of the star causing the lensing effect."
)
# Slider for the radius of the lensing star (Star 2)
radius_star2_solar = st.sidebar.slider(
    "Radius of Lensing Star (Star 2) (Solar Radii)",
    min_value=0.1, max_value=5.0, value=1.0, step=0.1,
    help="Adjust the visual size of the lensing star. (Does not affect lensing physics in this simplified model)."
)

# Slider for the distance of Star 2 from Earth (as a fraction of Earth-Star1 distance)
distance_lens_fraction = st.sidebar.slider(
    "Distance of Star 2 from Earth (Fraction of Earth-Star1 distance)",
    min_value=0.1, max_value=0.9, value=0.5, step=0.05,
    help="Position Star 2 between Earth and Star 1. 0.5 means halfway."
)

# Checkbox to decide if Star 2 has an orbiting planet
has_planet = st.sidebar.checkbox("Star 2 has an orbiting Planet?", value=False)

mass_planet_jupiter = 0.0
radius_planet_jupiter = 0.0
if has_planet:
    # Sliders for planet mass and radius, visible only if has_planet is checked
    mass_planet_jupiter = st.sidebar.slider(
        "Mass of Planet (Jupiter Masses)",
        min_value=0.1, max_value=5.0, value=1.0, step=0.1,
        help="Adjust the mass of the orbiting planet."
    )
    radius_planet_jupiter = st.sidebar.slider(
        "Radius of Planet (Jupiter Radii)",
        min_value=0.1, max_value=2.0, value=1.0, step=0.1,
        help="Adjust the visual size of the orbiting planet."
    )
    st.sidebar.info("Note: Planet's orbit is fixed for simplicity. Its lensing effect is a simplified perturbation.")

st.sidebar.write("---")
start_simulation = st.sidebar.button("🚀 Start Simulation")

# --- Simulation Setup ---

# Define fixed conceptual distances for the simulation visualization
# These are not real-world units but represent relative distances for plotting
D_S = 1000  # Conceptual distance from Earth to Star 1 (Source)
D_L = D_S * distance_lens_fraction  # Conceptual distance from Earth to Star 2 (Lens)
D_LS = D_S - D_L # Conceptual distance from Star 2 to Star 1

# Convert slider values to SI units for physical calculations
M_lens_SI = mass_star2_solar * M_SUN
R_lens_SI = radius_star2_solar * R_SUN # Corrected to use R_SUN
M_planet_SI = mass_planet_jupiter * M_JUPITER
R_planet_SI = radius_planet_jupiter * R_JUPITER # Corrected to use R_JUPITER

# --- Gravitational Lensing Physics Functions ---

def einstein_radius(M, D_L_conceptual, D_S_conceptual, D_LS_conceptual):
    """
    Calculates the Einstein radius in radians.
    Uses conceptual distances for the ratio, but SI mass.
    """
    if D_L_conceptual == 0 or D_S_conceptual == 0 or D_LS_conceptual == 0:
        return 0.0
    
    scaling_factor_for_vis = 1e10 
    
    return np.sqrt((4 * G * M / c**2) * (D_LS_conceptual / (D_L_conceptual * D_S_conceptual))) * scaling_factor_for_vis

def magnification(u):
    """
    Calculates the magnification for a point source/lens.
    u is the dimensionless impact parameter (beta / Einstein_radius).
    """
    if u == 0: 
        return 100.0 
    return (u**2 + 2) / (u * np.sqrt(u**2 + 4))

# --- Simulation Visualization Area ---
col1, col2 = st.columns([2, 1]) 

with col1:
    st.header("Visual Simulation")
    simulation_placeholder = st.empty()

with col2:
    st.header("Brightness Curve (Magnification)")
    brightness_placeholder = st.empty()

# --- Run Simulation on Button Click ---
if start_simulation:
    st.toast("Simulation started! ✨", icon="�")

    # Simulation parameters for the animation loop
    total_steps = 200 
    time_duration = 20 # Increased duration for more stability
    dt = time_duration / total_steps 

    cross_distance_max = 0.05 * D_S

    times = []
    magnifications_data = []

    fig_sim, ax_sim = plt.subplots(figsize=(10, 5))
    fig_bright, ax_bright = plt.subplots(figsize=(8, 4))

    # Configure simulation plot (ax_sim) for a space theme
    ax_sim.set_facecolor('black') 
    fig_sim.patch.set_facecolor('black') 
    ax_sim.set_xlim(-0.1 * D_S, D_S * 1.1) 
    ax_sim.set_ylim(-cross_distance_max * 1.5, cross_distance_max * 1.5) 
    ax_sim.set_aspect('equal', adjustable='box') 
    ax_sim.axis('off') 

    # Configure brightness plot (ax_bright)
    ax_bright.set_facecolor('black') 
    fig_bright.patch.set_facecolor('black') 
    ax_bright.set_xlabel("Time (conceptual units)", color='white') 
    ax_bright.set_ylabel("Magnification", color='white') 
    ax_bright.tick_params(axis='x', colors='white') 
    ax_bright.tick_params(axis='y', colors='white') 
    ax_bright.spines['bottom'].set_color('white')
    ax_bright.spines['left'].set_color('white')
    ax_bright.spines['top'].set_color('white')
    ax_bright.spines['right'].set_color('white')
    ax_bright.set_title("Apparent Brightness of Star 1", color='white') 
    ax_bright.grid(True, linestyle='--', alpha=0.5, color='gray') 

    # --- Animation Loop ---
    for i in range(total_steps):
        current_time = i * dt
        y_lens = -cross_distance_max + (2 * cross_distance_max / total_steps) * i
        x_lens = D_L 

        x_planet, y_planet = x_lens, y_lens 
        if has_planet:
            planet_orbit_radius = 0.015 * D_S 
            planet_angle = (current_time / time_duration) * 2 * np.pi * 5
            x_planet = x_lens + planet_orbit_radius * np.cos(planet_angle)
            y_planet = y_lens + planet_orbit_radius * np.sin(planet_angle)

        beta_star2 = np.abs(y_lens) / D_L 
        theta_E_star2 = einstein_radius(M_lens_SI, D_L, D_S, D_LS)
        u_star2 = beta_star2 / theta_E_star2 if theta_E_star2 > 0 else float('inf')
        mag_star2 = magnification(u_star2)

        mag_planet_effect = 1.0 
        if has_planet:
            beta_planet = np.abs(y_planet) / D_L 
            theta_E_planet = einstein_radius(M_planet_SI, D_L, D_S, D_LS)
            u_planet = beta_planet / theta_E_planet if theta_E_planet > 0 else float('inf')
            
            if u_planet < 2:
                mag_planet_effect = magnification(u_planet)
            else:
                mag_planet_effect = 1.0 

        total_magnification = mag_star2 * mag_planet_effect
        
        times.append(current_time)
        magnifications_data.append(total_magnification)

        # --- Update Visual Simulation Plot ---
        ax_sim.clear() 
        ax_sim.set_facecolor('black') 
        ax_sim.set_xlim(-0.1 * D_S, D_S * 1.1)
        ax_sim.set_ylim(-cross_distance_max * 1.5, cross_distance_max * 1.5)
        ax_sim.set_aspect('equal', adjustable='box')
        ax_sim.axis('off')

        ax_sim.plot(0, 0, 'o', color='blue', markersize=10, label="Earth (Observer)")
        ax_sim.text(0, -0.08 * D_S, "Earth", color='white', ha='center', fontsize=10)

        ax_sim.plot(D_S, 0, '*', color='yellow', markersize=15, label="Star 1 (Source)")
        ax_sim.text(D_S, -0.08 * D_S, "Star 1", color='white', ha='center', fontsize=10)

        star2_circle = plt.Circle((x_lens, y_lens), radius=R_lens_SI * 2e-9, color='red', alpha=0.8)
        ax_sim.add_patch(star2_circle)
        ax_sim.text(x_lens, y_lens + R_lens_SI * 2e-9 + 0.005 * D_S, "Star 2", color='white', ha='center', fontsize=10)

        if has_planet:
            planet_circle = plt.Circle((x_planet, y_planet), radius=R_planet_SI * 2e-9, color='cyan', alpha=0.8)
            ax_sim.add_patch(planet_circle)
            ax_sim.text(x_planet, y_planet + R_planet_SI * 2e-9 + 0.005 * D_S, "Planet", color='white', ha='center', fontsize=10)

        if u_star2 < 5: 
            ax_sim.plot([D_S, x_lens], [0, y_lens], 'w--', alpha=0.5) 
            ax_sim.plot([x_lens, 0], [y_lens, 0], 'w--', alpha=0.5) 
            ax_sim.arrow(D_S, 0, x_lens - D_S, y_lens, head_width=0.005*D_S, head_length=0.01*D_S, fc='white', ec='white', length_includes_head=True)
            ax_sim.arrow(x_lens, y_lens, 0 - x_lens, 0 - y_lens, head_width=0.005*D_S, head_length=0.01*D_S, fc='white', ec='white', length_includes_head=True)
        else:
            ax_sim.plot([D_S, 0], [0, 0], 'w-', alpha=0.5) 
            ax_sim.arrow(D_S, 0, -D_S, 0, head_width=0.005*D_S, head_length=0.01*D_S, fc='white', ec='white', length_includes_head=True)

        simulation_placeholder.pyplot(fig_sim)

        # --- Update Brightness Curve Plot ---
        ax_bright.clear() 
        ax_bright.set_facecolor('black') 
        ax_bright.plot(times, magnifications_data, color='lime', linewidth=2) 
        ax_bright.set_xlabel("Time (conceptual units)", color='white')
        ax_bright.set_ylabel("Magnification", color='white')
        ax_bright.tick_params(axis='x', colors='white')
        ax_bright.tick_params(axis='y', colors='white')
        ax_bright.spines['bottom'].set_color('white')
        ax_bright.spines['left'].set_color('white')
        ax_bright.spines['top'].set_color('white')
        ax_bright.spines['right'].set_color('white')
        ax_bright.set_title("Apparent Brightness of Star 1", color='white')
        ax_bright.grid(True, linestyle='--', alpha=0.5, color='gray')
        ax_bright.set_ylim(0.8, max(magnifications_data) * 1.2 if magnifications_data else 2.0)

        brightness_placeholder.pyplot(fig_bright)

        time.sleep(0.08) # Increased delay for more stability

    plt.close(fig_sim)
    plt.close(fig_bright)

    st.toast("Simulation finished! 🎉", icon="✅")

# --- Footer / Instructions ---
st.markdown(
    """
    <div style="text-align: center; margin-top: 30px; color: gray;">
        <p>This simulator provides a simplified visualization of gravitational lensing.</p>
        <p>Adjust the parameters in the sidebar to observe how the mass and position of a foreground star (Star 2) and its planet affect the apparent brightness of a background star (Star 1).</p>
    </div>
    """,
    unsafe_allow_html=True
)
