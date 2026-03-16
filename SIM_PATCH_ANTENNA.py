# -*- coding: utf-8 -*-
import os, glob
import numpy as np
import matplotlib.pyplot as plt

# 1. PARAMETER MANAGEMENT

# Free Parameters (Geometric Sweep)
f0_list = [1.7475e9]            # Target center frequency (Hz)
D_list = [34.5]                 # Length of each patch side (mm)
L_cut_list = [3.0]              # Corner cut length for circular polarization (mm)
divisors_list = [2.0]           # Wavelength divisor for patch spacing
thickness_list = [1.6]          # FR4 substrate height (mm)
alpha_deg_list = [90.0]         # Angle between cross arms (90 = orthogonal cross)

# Multi-Fidelity Parameters
t_copper = 0.0                  # Copper thickness (mm). 0.0 = Fast (2D). 0.035 = Precise/Manufacturing (3D, 1 oz).
max_timesteps = 60000           # FDTD Limit. Increase to 150000 if t_copper > 0 to avoid premature cuts.

# Fixed Parameters (Materials and Internal Geometry)
epsR_substrate = 4.5            # FR4 relative permittivity
epsR_connector = 2.1            # Teflon permittivity for SMA connector
feed_resistance = 50.0          # Port impedance (Ohms)
feed_offset_x = 7.0             # Distance from center inwards for 50 ohms matching (mm)
outer_pin_length = 5.0          # SMA connector pin length (mm)
sma_outer_radius = 2.5          # SMA shield outer radius (mm)
sma_inner_radius = 2.0          # SMA inner pin radius (mm)

# Fixed Parameters (FDTD Simulation and Mesh)
c0 = 299792458
edge_margin = 40.0              # Extra substrate margin from the furthest antenna (mm)
air_box_margin = 45.0           # Air above the antenna before the PML limit (mm)
substrate_cells = 4             # Vertical discretization in the substrate
beam_drop_db = 3.0              # dB drop from maximum to measure physical lobe width

# Flow Control Flags
post_process_only = False       # Flag to avoid simulating if data already exists
verify_first_geometry = True    # Flag to open the 3D interface before starting the sweep


def calculate_hpbw(theta, pattern_dB, drop_db): 
    max_idx = np.argmax(pattern_dB)
    max_val = pattern_dB[max_idx]
    limit_db = max_val - drop_db
    
    left = np.where(pattern_dB[:max_idx] <= limit_db)[0]
    theta_left = theta[left[-1]] if len(left) > 0 else theta[0]
    right = np.where(pattern_dB[max_idx:] <= limit_db)[0]
    theta_right = theta[max_idx + right[0]] if len(right) > 0 else theta[-1]
    
    return abs(theta_right - theta_left)

def calculate_rmse(theta, spectrum_dB, target_angles):
    """ Calculates the Root Mean Square Error of the peak localization """
    error_sq = 0.0
    for ang in target_angles:
        range_idx = np.abs(theta - ang) < 20
        if not np.any(range_idx): return 999.0
        peak_idx = np.argmax(spectrum_dB[range_idx])
        calc_ang = theta[range_idx][peak_idx]
        error_sq += (calc_ang - ang)**2
    return np.sqrt(error_sq / len(target_angles))

openems_path = r'C:\Users\herce\openEMS\openEMS'
if hasattr(os, 'add_dll_directory'): os.add_dll_directory(openems_path)
os.environ['PATH'] = openems_path + os.pathsep + os.environ.get('PATH', '')

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

final_results = []

# 2. SAVING CONFIGURATION
root_output = r'C:\Users\herce\Conda\Sims\Simple_Patch_Antenna' 
base_sweep = os.path.join(root_output, "Auto_Sweep")

sweep_counter = 1
sweep_folder = f"{base_sweep}_{sweep_counter}"
while os.path.exists(sweep_folder):
    sweep_counter += 1
    sweep_folder = f"{base_sweep}_{sweep_counter}"

os.makedirs(sweep_folder)
report_file = os.path.join(sweep_folder, "Sweep_Report.txt")

with open(report_file, "w") as f:
    f.write(f"{'f0 (GHz)':<9} | {'D (mm)':<7} | {'L (mm)':<7} | {'Div':<5} | {'h (mm)':<7} | {'f_res (GHz)':<12} | {'S11 (dB)':<9} | {'SWR':<5} | {'CT (dB)':<8} | {'FoV_E(°)':<9} | {'Res_B(°)':<9} | {'Res_M(°)':<9} | {'RMSE(°)':<9}\n")
    f.write("-" * 138 + "\n")

# DSP MODULE: DIRECTION OF ARRIVAL (MUSIC AND BARTLETT)

def simulate_and_process_doa(manifold, theta_array, source_angles, SNR_dB=15):
    """
    Generates a standardized DOA scenario and calculates pseudo-spectrums.
    manifold: Complex matrix (N_antennas, N_angles) exported from openEMS.
    """
    K = len(source_angles)
    M, N_theta = manifold.shape
    
    # 1. Build covariance matrix
    A_sources = []
    for ang in source_angles:
        idx = np.argmin(np.abs(theta_array - ang))
        A_sources.append(manifold[:, idx])
    A_sources = np.column_stack(A_sources)
    
    np.random.seed(42) # Fixed seed
    snapshots = 200
    S = (np.random.randn(K, snapshots) + 1j * np.random.randn(K, snapshots)) / np.sqrt(2)
    S = S * (10**(SNR_dB/20))
    X = A_sources @ S
    
    noise = (np.random.randn(M, snapshots) + 1j * np.random.randn(M, snapshots)) / np.sqrt(2)
    X = X + noise
    Rx = (X @ X.conj().T) / snapshots
    
    # 2. MUSIC Algorithm
    vals, vecs = np.linalg.eigh(Rx)
    En = vecs[:, :-K] # Noise subspace
    
    P_music = np.zeros(N_theta)
    P_bartlett = np.zeros(N_theta)
    
    for i in range(N_theta):
        a = manifold[:, i]
        # MUSIC
        den = np.abs(a.conj().T @ En @ En.conj().T @ a)
        P_music[i] = 1.0 / (den + 1e-12)
        # Bartlett
        P_bartlett[i] = np.abs(a.conj().T @ Rx @ a)
        
    P_music = 10 * np.log10(P_music / np.max(P_music))
    P_bartlett = 10 * np.log10(P_bartlett / np.max(P_bartlett))
    
    return P_music, P_bartlett

def measure_doa_peak_width(theta, spectrum_dB, expected_angle, drop=3.0):
    """ Measures peak sharpness (Resolution Power). Returns 999 if unresolved. """
    range_idx = np.where((theta >= expected_angle - 10) & (theta <= expected_angle + 10))[0]
    if len(range_idx) == 0: return 999.0
    
    local_peak_idx = range_idx[np.argmax(spectrum_dB[range_idx])]
    max_val = spectrum_dB[local_peak_idx]
    limit = max_val - drop
    
    left = local_peak_idx
    while left > 0 and spectrum_dB[left] > limit: left -= 1
    right = local_peak_idx
    while right < len(theta)-1 and spectrum_dB[right] > limit: right += 1
    
    width = theta[right] - theta[left]
    if width > 35: return 999.0
    return width



# 3. MASTER LOOP
first_simulation = True       

for f0 in f0_list:
    lambda_0 = (c0 / f0) * 1000
    for D in D_list:
        for L_cut in L_cut_list:
            for divisor in divisors_list:
                for thickness in thickness_list:
                    for alpha_deg in alpha_deg_list:
                        
                        element_spacing = lambda_0 / divisor
                        study_folder = f"{int(f0/1e6)}MHz_L{L_cut}_D{D}_Div{divisor}_A{int(alpha_deg)}"
                        study_path = os.path.join(sweep_folder, study_folder)
                        if not os.path.exists(study_path): os.makedirs(study_path)

                        n_sim = len(glob.glob(os.path.join(study_path, 'sim_*'))) + 1
                        Sim_Path = os.path.join(study_path, f"sim_{n_sim}")
                        os.makedirs(Sim_Path)

                        print(f"\n>>> STARTING FLOW: D={D}, L={L_cut}, Alpha={alpha_deg}°")
                        
                        alpha_rad = np.deg2rad(alpha_deg)
                        positions = [
                            (0, 0),                                                                         
                            (element_spacing, 0),                                                      
                            (-element_spacing, 0),                                                     
                            (element_spacing * np.cos(alpha_rad), element_spacing * np.sin(alpha_rad)),    
                            (-element_spacing * np.cos(alpha_rad), -element_spacing * np.sin(alpha_rad))   
                        ]

                        # Residual floating point cleanup to prevent mesh errors
                        positions = [(np.round(x, 4), np.round(y, 4)) for x, y in positions]

                        max_x = max([abs(x) for x, y in positions]) + (D / 2.0)
                        max_y = max([abs(y) for x, y in positions]) + (D / 2.0)
                        substrate_width = (2.0 * max_x) + edge_margin
                        substrate_length = (2.0 * max_y) + edge_margin
                        kappa = 1e-3 * 2*np.pi*2.45e9 * EPS0 * epsR_substrate  

                        # Global Variables for this Parameter Set
                        theta_doa = np.arange(-180.0, 180.0, 1.0)
                        manifold_complex = np.zeros((5, len(theta_doa)), dtype=complex)
                        
                        calc_f_res = f0
                        min_s11 = 0.0
                        min_swr = 0.0
                        max_cross_talk = 0.0
                        hpbw_E, hpbw_H = 0.0, 0.0

                        # 4. SEQUENTIAL EXCITATION LOOP (x5)

                        for active_port in range(5):
                            Sim_Path_Port = os.path.join(Sim_Path, f'port_{active_port}')
                            os.makedirs(Sim_Path_Port)
                            
                            print(f"    -> Simulating FDTD Active Port: {active_port}/4 ...")

                            SimBox = np.array([substrate_width + 2*air_box_margin, substrate_length + 2*air_box_margin, 200])
                            fc = 1e9                                

                            FDTD = openEMS(NrTS=max_timesteps, EndCriteria=1e-4) # Dynamic limit
                            FDTD.SetGaussExcite(f0, fc)
                            FDTD.SetBoundaryCond(['PML_8']*6)

                            CSX = ContinuousStructure()
                            FDTD.SetCSX(CSX)
                            mesh = CSX.GetGrid()
                            mesh.SetDeltaUnit(1e-3)
                            mesh_res = C0/(f0+fc)/1e-3/20   

                            mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
                            mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
                            
                            # Conditional Z mesh for copper thickness
                            z_lines = np.linspace(0, thickness, substrate_cells+1).tolist()
                            if t_copper > 0.0:
                                z_lines.append(thickness + t_copper)
                            mesh.AddLine('z', z_lines)
                            mesh.AddLine('z', [-SimBox[2]/3, SimBox[2]*2/3])

                            patch = CSX.AddMetal('patch')
                            ports = []

                            for i, (x_pos, y_pos) in enumerate(positions):
                                points_x = [x_pos - D/2, x_pos - D/2, x_pos - D/2 + L_cut, x_pos + D/2, x_pos + D/2, x_pos + D/2 - L_cut]
                                points_y = [y_pos + D/2, y_pos - D/2 + L_cut, y_pos - D/2, y_pos - D/2, y_pos + D/2 - L_cut, y_pos + D/2]
                                
                                # Conditional 3D Extrusion
                                if t_copper == 0.0:
                                    patch.AddPolygon([points_x, points_y], 'z', thickness, priority=10)
                                else:
                                    patch.AddLinPoly([points_x, points_y], 'z', thickness, t_copper, priority=10)
                                
                                sma_dielectric = CSX.AddMaterial(f'sma_dielectric_{i}', epsilon=epsR_connector)
                                sma_shield = CSX.AddMetal(f'sma_shield_{i}')
                                
                                sma_shield.AddBox([x_pos - feed_offset_x - sma_outer_radius, y_pos - sma_outer_radius, -outer_pin_length], 
                                                  [x_pos - feed_offset_x + sma_outer_radius, y_pos + sma_outer_radius, 0], priority=1)
                                sma_dielectric.AddBox([x_pos - feed_offset_x - sma_inner_radius, y_pos - sma_inner_radius, -outer_pin_length], 
                                                       [x_pos - feed_offset_x + sma_inner_radius, y_pos + sma_inner_radius, 0], priority=2)

                                port_start = [x_pos - feed_offset_x, y_pos, 0]
                                port_end  = [x_pos - feed_offset_x, y_pos, thickness]
                                
                                # Only injects 1V to the active port, others 0V (Passive)
                                excitation_voltage = 1.0 if i == active_port else 0.0
                                p = FDTD.AddLumpedPort(i+1, feed_resistance, port_start, port_end, 'z', excitation_voltage, priority=5, edges2grid='xy')
                                ports.append(p)

                            substrate = CSX.AddMaterial('substrate', epsilon=epsR_substrate, kappa=kappa)
                            substrate.AddBox([-substrate_width/2, -substrate_length/2, 0], [substrate_width/2, substrate_length/2, thickness], priority=0)

                            ground_plane = CSX.AddMetal('ground_plane')
                            ground_plane.AddBox([-substrate_width/2, -substrate_length/2, 0], [substrate_width/2, substrate_length/2, 0], priority=10)

                            FDTD.AddEdges2Grid(dirs='xy', properties=patch, metal_edge_res=mesh_res/2)
                            mesh.SmoothMeshLines('all', mesh_res, 1.4)

                            nf2ff = FDTD.CreateNF2FFBox()
                            
                            csx_file = os.path.join(Sim_Path_Port, 'simp_patch_array.xml')
                            CSX.Write2XML(csx_file)

                            if active_port == 0 and verify_first_geometry and first_simulation:
                                print("\n>>> OPENING AppCSXCAD FOR VISUAL VERIFICATION...")
                                os.system(f'AppCSXCAD "{csx_file}"')
                                response = input("Continue auto sweep? (y/n): ")
                                if response.lower() != 'y': import sys; sys.exit()
                                first_simulation = False

                            if not post_process_only:
                                FDTD.Run(Sim_Path_Port, verbose=0, cleanup=True)

                            # 5. POST-PROCESS (HARDWARE PHYSICS)

                            if active_port == 0:
                                f = np.linspace(max(1e9, f0-fc), f0+fc, 401)
                                ports[0].CalcPort(Sim_Path_Port, f)
                                
                                s11_dB = 20.0 * np.log10(np.abs(ports[0].uf_ref / (ports[0].uf_inc + 1e-15)))
                                f_res_idx = np.argmin(s11_dB)
                                min_s11 = s11_dB[f_res_idx]
                                calc_f_res = f[f_res_idx] # FIX RESONANT FREQUENCY FOR ALL

                                # SWR Calculation
                                gamma = 10**(min_s11 / 20.0)
                                min_swr = (1 + np.abs(gamma)) / (1 - np.abs(gamma))

                                # Cross Talk (Mutual Coupling) Calculation
                                couplings = []
                                for p_idx in range(1, 5):
                                    ports[p_idx].CalcPort(Sim_Path_Port, f)
                                    s_ix_dB = 20.0 * np.log10(np.abs(ports[p_idx].uf_ref / (ports[0].uf_inc + 1e-15)))
                                    couplings.append(s_ix_dB[f_res_idx]) # Extracts coupling at exact f_res
                                max_cross_talk = np.max(couplings) # Worst case coupling

                                plt.figure(figsize=(8,4))
                                plt.plot(f/1e9, s11_dB, 'k-')
                                plt.grid(True); plt.ylabel('S11 (dB)'); plt.xlabel('Frequency (GHz)')
                                plt.savefig(os.path.join(Sim_Path, '1_S11_Central.png')); plt.close()

                                # Far-Field of Port 0 (Hardware FoV)
                                nf2ff_res_0 = nf2ff.CalcNF2FF(Sim_Path_Port, calc_f_res, theta_doa, [0., 90.], center=[0,0,0])
                                fields_matrix = nf2ff_res_0.E_norm[0]
                                
                                E_plane_norm = np.squeeze(20.0*np.log10(fields_matrix[:, 0] / (np.max(fields_matrix[:, 0]) + 1e-15)))
                                H_plane_norm = np.squeeze(20.0*np.log10(fields_matrix[:, 1] / (np.max(fields_matrix[:, 1]) + 1e-15)))
                                
                                hpbw_E = calculate_hpbw(theta_doa, E_plane_norm, beam_drop_db)
                                hpbw_H = calculate_hpbw(theta_doa, H_plane_norm, beam_drop_db)

                                plt.figure(figsize=(8,4))
                                plt.plot(theta_doa, E_plane_norm + nf2ff_res_0.Dmax[0], 'k-', label=f'E-plane FoV: {hpbw_E:.1f}°')
                                plt.plot(theta_doa, H_plane_norm + nf2ff_res_0.Dmax[0], 'r--', label=f'H-plane FoV: {hpbw_H:.1f}°')
                                plt.axhline(nf2ff_res_0.Dmax[0] - beam_drop_db, color='gray', linestyle=':')
                                plt.grid(True); plt.legend(); plt.ylabel('Directivity (dBi)'); plt.xlabel('Theta Angle (Degrees)')
                                plt.savefig(os.path.join(Sim_Path, '2_FoV_Fisico.png'))
                                plt.close() # Closes figure so it doesn't print to console

                            # Strict Complex Phase Extraction referenced to [0,0,0]
                            nf2ff_calib = nf2ff.CalcNF2FF(Sim_Path_Port, calc_f_res, theta_doa, [0.], center=[0,0,0])
                            manifold_complex[active_port, :] = nf2ff_calib.E_theta[0][0, :] # Co-polar complex component

                        # 6. DIGITAL POST-PROCESS (DOA RESOLUTION)

                        print("    -> Calculating (MUSIC / Bartlett) ...")
                        
                        # MANIFOLD NORMALIZATION (For SNR and synthetic noise to work properly)
                        manifold_complex = manifold_complex / np.max(np.abs(manifold_complex))
                        
                        escenario_angles = [-25.0, 25.0]
                        P_music, P_bartlett = simulate_and_process_doa(manifold_complex, theta_doa, escenario_angles, SNR_dB=15)
                        
                        res_bartlett = measure_doa_peak_width(theta_doa, P_bartlett, escenario_angles[1])
                        res_music = measure_doa_peak_width(theta_doa, P_music, escenario_angles[1])
                        rmse_music = calculate_rmse(theta_doa, P_music, escenario_angles)

                        # Software spectrum plot
                        plt.figure(figsize=(10,5))
                        plt.plot(theta_doa, P_bartlett, 'b--', label=f'Bartlett (Res: {res_bartlett:.1f}°)')
                        plt.plot(theta_doa, P_music, 'r-', lw=2, label=f'MUSIC (Res: {res_music:.1f}°)')
                        plt.axvline(-10, color='k', linestyle=':', alpha=0.5)
                        plt.axvline(10, color='k', linestyle=':', alpha=0.5)
                        plt.title(f'DOA Resolution Test (Zenith to $\pm10^\circ$)')
                        plt.xlabel('Theta Angle (Degrees)'); plt.ylabel('Pseudo-Spectrum (dB)')
                        plt.grid(True); plt.legend(); plt.ylim([-40, 5])
                        plt.savefig(os.path.join(Sim_Path, '3_Espectro_DOA.png')); plt.close()

                        # 7. FINAL REPORT

                        res_B_str = f"{res_bartlett:.1f}" if res_bartlett != 999.0 else "FAIL"
                        res_M_str = f"{res_music:.1f}" if res_music != 999.0 else "FAIL"
                        rmse_M_str = f"{rmse_music:.2f}" if rmse_music != 999.0 else "FAIL"

                        config = {
                            'f0_GHz': float(f0/1e9), 'D_mm': float(D), 'L_mm': float(L_cut),
                            'Divisor': float(divisor), 'h_mm': float(thickness),
                            'f_res_GHz': float(round(calc_f_res/1e9, 3)), 'S11_dB': float(round(min_s11, 2)),
                            'SWR': float(round(min_swr, 2)), 'CT_dB': float(round(max_cross_talk, 2)),
                            'FoV_E': float(round(hpbw_E, 1)), 'Res_B': res_B_str, 'Res_M': res_M_str, 'RMSE': rmse_M_str
                        }
                        final_results.append(config)
                        
                        with open(report_file, "a") as f:
                            f.write(f"{config['f0_GHz']:<9.3f} | {config['D_mm']:<7.1f} | {config['L_mm']:<7.1f} | {config['Divisor']:<5.1f} | {config['h_mm']:<7.1f} | {config['f_res_GHz']:<12.3f} | {config['S11_dB']:<9.2f} | {config['SWR']:<5.2f} | {config['CT_dB']:<8.2f} | {config['FoV_E']:<9.1f} | {config['Res_B']:<9} | {config['Res_M']:<9} | {config['RMSE']:<9}\n")

# 8. FINAL CONSOLE PRINT
print("\n" + "="*138)
print(f"{'f0 (GHz)':<9} | {'D (mm)':<7} | {'L (mm)':<7} | {'Div':<5} | {'h (mm)':<7} | {'f_res (GHz)':<12} | {'S11 (dB)':<9} | {'SWR':<5} | {'CT (dB)':<8} | {'FoV_E(°)':<9} | {'Res_B(°)':<9} | {'Res_M(°)':<9} | {'RMSE(°)':<9}")
print("-" * 138)
for r in final_results:
    print(f"{r['f0_GHz']:<9.3f} | {r['D_mm']:<7.1f} | {r['L_mm']:<7.1f} | {r['Divisor']:<5.1f} | {r['h_mm']:<7.1f} | {r['f_res_GHz']:<12.3f} | {r['S11_dB']:<9.2f} | {r['SWR']:<5.2f} | {r['CT_dB']:<8.2f} | {r['FoV_E']:<9.1f} | {r['Res_B']:<9} | {r['Res_M']:<9} | {r['RMSE']:<9}")
print("="*138 + "\n")
print(f"Everything saved in: {sweep_folder}")