# -*- coding: utf-8 -*-
import os, glob
import numpy as np
import matplotlib.pyplot as plt

# 1. GESTIÓN DE PARÁMETROS

# Parámetros Libres (Barrido Geométrico)
lista_f0 = [1.7475e9]           # Frecuencia central(Hz)
lista_D = [34.5]                # Largo de cada lado del parche (mm)
lista_L_corte = [3.0]           # Largo del corte en las esquinas para polarización circular (mm)
lista_divisores = [2.0]         # Divisor de longitud de onda para separación entre parches
lista_espesores = [1.6]         # Altura del sustrato de FR4 (mm)
lista_alpha_deg = [90.0]        # Ángulo entre brazos de la cruz (90 = cruz ortogonal)

# Parámetros de Multi-Fidelidad
t_cobre = 0.0                   # Espesor del cobre (mm). 0.0 = Rápido (2D). 0.035 = Preciso/Fabricación (3D, 1 oz).
max_timesteps = 60000           # Límite FDTD. Subir a 150000 si t_cobre > 0 para evitar cortes prematuros.

# Parámetros Fijos (Materiales y Geometría Interna)
epsR_sustrato = 4.5             # Permitividad relativa del FR4
epsR_conector = 2.1             # Permitividad del teflón para el conector SMA
resistencia_alimentacion = 50.0 # Impedancia del puerto (Ohms)
desplazamiento_feed_x = 7.0     # Distancia desde el centro hacia adentro para adaptar a 50 ohms (mm)
largo_pin_exterior = 5.0        # Largo del pin del conector SMA (mm)
radio_exterior_sma = 2.5        # Radio exterior del escudo del SMA (mm)
radio_interior_sma = 2.0        # Radio del pin interior del SMA (mm)

# Parámetros Fijos (Simulación FDTD y Malla)
c0 = 299792458
margen_borde = 40.0             # Borde extra de sustrato desde la antena más lejana (mm)
colchon_aire = 45.0             # Aire sobre la antena antes del límite PML (mm)
celdas_sustrato = 4             # Discretización vertical en el sustrato
caida_haz_db = 3.0              # A cuántos dB del máximo medimos el ancho del lóbulo físico

# Flags de Control de Flujo
solo_post_proceso = False           # Flag para evitar simular si ya existen los datos
verificar_primera_geometria = True  # Flag para abrir la interfaz 3D antes de arrancar el barrido


def calcular_hpbw(theta, patron_dB, caida_db): 
    max_idx = np.argmax(patron_dB)
    max_val = patron_dB[max_idx]
    limite_db = max_val - caida_db
    
    izq = np.where(patron_dB[:max_idx] <= limite_db)[0]
    theta_izq = theta[izq[-1]] if len(izq) > 0 else theta[0]
    der = np.where(patron_dB[max_idx:] <= limite_db)[0]
    theta_der = theta[max_idx + der[0]] if len(der) > 0 else theta[-1]
    
    return abs(theta_der - theta_izq)

def calcular_rmse(theta, espectro_dB, angulos_reales):
    """ Calcula el Error Cuadrático Medio de la localización de los picos """
    error_sq = 0.0
    for ang in angulos_reales:
        rango = np.abs(theta - ang) < 20
        if not np.any(rango): return 999.0
        idx_pico = np.argmax(espectro_dB[rango])
        ang_calc = theta[rango][idx_pico]
        error_sq += (ang_calc - ang)**2
    return np.sqrt(error_sq / len(angulos_reales))

openems_path = r'C:\Users\herce\openEMS\openEMS'
if hasattr(os, 'add_dll_directory'): os.add_dll_directory(openems_path)
os.environ['PATH'] = openems_path + os.pathsep + os.environ.get('PATH', '')

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

resultados_finales = []

# 2. CONFIGURACIÓN DE GUARDADO
root_output = r'C:\Users\herce\Conda\Sims\Simple_Patch_Antenna' 
base_barrido = os.path.join(root_output, "Barrido_Automatico")

contador_barrido = 1
carpeta_barrido = f"{base_barrido}_{contador_barrido}"
while os.path.exists(carpeta_barrido):
    contador_barrido += 1
    carpeta_barrido = f"{base_barrido}_{contador_barrido}"

os.makedirs(carpeta_barrido)
archivo_reporte = os.path.join(carpeta_barrido, "Reporte_Barrido.txt")

with open(archivo_reporte, "w") as f:
    f.write(f"{'f0 (GHz)':<9} | {'D (mm)':<7} | {'L (mm)':<7} | {'Div':<5} | {'h (mm)':<7} | {'f_res (GHz)':<12} | {'S11 (dB)':<9} | {'SWR':<5} | {'CT (dB)':<8} | {'FoV_E(°)':<9} | {'Res_B(°)':<9} | {'Res_M(°)':<9} | {'RMSE(°)':<9}\n")
    f.write("-" * 138 + "\n")

# MÓDULO DSP: DIRECCIÓN DE ARRIBO (MUSIC Y BARTLETT)

def simular_y_procesar_doa(manifold, theta_array, angulos_fuentes, SNR_dB=15):
    """
    Genera un escenario DOA estandarizado y calcula los pseudo-espectros.
    manifold: Matriz compleja (N_antenas, N_angulos) exportada de openEMS.
    """
    K = len(angulos_fuentes)
    M, N_theta = manifold.shape
    
    # 1. Armar matriz de covarianza
    A_sources = []
    for ang in angulos_fuentes:
        idx = np.argmin(np.abs(theta_array - ang))
        A_sources.append(manifold[:, idx])
    A_sources = np.column_stack(A_sources)
    
    np.random.seed(42) # Dejar semilla fija
    snapshots = 200
    S = (np.random.randn(K, snapshots) + 1j * np.random.randn(K, snapshots)) / np.sqrt(2)
    S = S * (10**(SNR_dB/20))
    X = A_sources @ S
    
    ruido = (np.random.randn(M, snapshots) + 1j * np.random.randn(M, snapshots)) / np.sqrt(2)
    X = X + ruido
    Rx = (X @ X.conj().T) / snapshots
    
    # 2. Algoritmo MUSIC
    vals, vecs = np.linalg.eigh(Rx)
    En = vecs[:, :-K] # Subespacio de ruido
    
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

def medir_ancho_pico_doa(theta, espectro_dB, angulo_esperado, caida=3.0):
    """ Mide la agudeza del pico (Poder de Resolución). Devuelve 999 si no resuelve. """
    rango = np.where((theta >= angulo_esperado - 10) & (theta <= angulo_esperado + 10))[0]
    if len(rango) == 0: return 999.0
    
    idx_pico_local = rango[np.argmax(espectro_dB[rango])]
    max_val = espectro_dB[idx_pico_local]
    limite = max_val - caida
    
    izq = idx_pico_local
    while izq > 0 and espectro_dB[izq] > limite: izq -= 1
    der = idx_pico_local
    while der < len(theta)-1 and espectro_dB[der] > limite: der += 1
    
    ancho = theta[der] - theta[izq]
    if ancho > 35: return 999.0
    return ancho



# 3. BUCLE MAESTRO
primera_simulacion = True       

for f0 in lista_f0:
    lambda_0 = (c0 / f0) * 1000
    for D in lista_D:
        for L_corte in lista_L_corte:
            for lambda_divisor in lista_divisores:
                for espesor_sustrato in lista_espesores:
                    for alpha_deg in lista_alpha_deg:
                        
                        separacion_elementos = lambda_0 / lambda_divisor
                        study_folder = f"{int(f0/1e6)}MHz_L{L_corte}_D{D}_Div{lambda_divisor}_A{int(alpha_deg)}"
                        study_path = os.path.join(carpeta_barrido, study_folder)
                        if not os.path.exists(study_path): os.makedirs(study_path)

                        n_sim = len(glob.glob(os.path.join(study_path, 'sim_*'))) + 1
                        Sim_Path = os.path.join(study_path, f"sim_{n_sim}")
                        os.makedirs(Sim_Path)

                        print(f"\n>>> INICIANDO FLUJO: D={D}, L={L_corte}, Alfa={alpha_deg}°")
                        
                        alpha_rad = np.deg2rad(alpha_deg)
                        posiciones = [
                            (0, 0),                                                                         
                            (separacion_elementos, 0),                                                      
                            (-separacion_elementos, 0),                                                     
                            (separacion_elementos * np.cos(alpha_rad), separacion_elementos * np.sin(alpha_rad)),    
                            (-separacion_elementos * np.cos(alpha_rad), -separacion_elementos * np.sin(alpha_rad))   
                        ]

                        # Limpieza de punto flotante residual para evitar errores de malla
                        posiciones = [(np.round(x, 4), np.round(y, 4)) for x, y in posiciones]

                        max_x = max([abs(x) for x, y in posiciones]) + (D / 2.0)
                        max_y = max([abs(y) for x, y in posiciones]) + (D / 2.0)
                        ancho_sustrato = (2.0 * max_x) + margen_borde
                        largo_sustrato = (2.0 * max_y) + margen_borde
                        kappa = 1e-3 * 2*np.pi*2.45e9 * EPS0 * epsR_sustrato  

                        # Variables Globales para este Set de Parámetros
                        theta_doa = np.arange(-180.0, 180.0, 1.0)
                        manifold_complex = np.zeros((5, len(theta_doa)), dtype=complex)
                        
                        f_res_calculada = f0
                        s11_minimo = 0.0
                        swr_minimo = 0.0
                        max_cross_talk = 0.0
                        hpbw_E, hpbw_H = 0.0, 0.0

                        # 4. BUCLE DE EXCITACIÓN SECUENCIAL (x5)

                        for puerto_activo in range(5):
                            Sim_Path_Puerto = os.path.join(Sim_Path, f'puerto_{puerto_activo}')
                            os.makedirs(Sim_Path_Puerto)
                            
                            print(f"    -> Simulando FDTD Puerto Activo: {puerto_activo}/4 ...")

                            SimBox = np.array([ancho_sustrato + 2*colchon_aire, largo_sustrato + 2*colchon_aire, 200])
                            fc = 1e9                                

                            FDTD = openEMS(NrTS=max_timesteps, EndCriteria=1e-4) # Límite dinámico
                            FDTD.SetGaussExcite(f0, fc)
                            FDTD.SetBoundaryCond(['PML_8']*6)

                            CSX = ContinuousStructure()
                            FDTD.SetCSX(CSX)
                            mesh = CSX.GetGrid()
                            mesh.SetDeltaUnit(1e-3)
                            mesh_res = C0/(f0+fc)/1e-3/20   

                            mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
                            mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
                            
                            # Malla en Z condicional para el espesor del cobre
                            z_lines = np.linspace(0, espesor_sustrato, celdas_sustrato+1).tolist()
                            if t_cobre > 0.0:
                                z_lines.append(espesor_sustrato + t_cobre)
                            mesh.AddLine('z', z_lines)
                            mesh.AddLine('z', [-SimBox[2]/3, SimBox[2]*2/3])

                            parche = CSX.AddMetal('parche')
                            puertos = []

                            for i, (x_pos, y_pos) in enumerate(posiciones):
                                puntos_x = [x_pos - D/2, x_pos - D/2, x_pos - D/2 + L_corte, x_pos + D/2, x_pos + D/2, x_pos + D/2 - L_corte]
                                puntos_y = [y_pos + D/2, y_pos - D/2 + L_corte, y_pos - D/2, y_pos - D/2, y_pos + D/2 - L_corte, y_pos + D/2]
                                
                                # Extrusión 3D condicional
                                if t_cobre == 0.0:
                                    parche.AddPolygon([puntos_x, puntos_y], 'z', espesor_sustrato, priority=10)
                                else:
                                    parche.AddLinPoly([puntos_x, puntos_y], 'z', espesor_sustrato, t_cobre, priority=10)
                                
                                dielectrico_sma = CSX.AddMaterial(f'dielectrico_sma_{i}', epsilon=epsR_conector)
                                escudo_sma = CSX.AddMetal(f'escudo_sma_{i}')
                                
                                escudo_sma.AddBox([x_pos - desplazamiento_feed_x - radio_exterior_sma, y_pos - radio_exterior_sma, -largo_pin_exterior], 
                                                  [x_pos - desplazamiento_feed_x + radio_exterior_sma, y_pos + radio_exterior_sma, 0], priority=1)
                                dielectrico_sma.AddBox([x_pos - desplazamiento_feed_x - radio_interior_sma, y_pos - radio_interior_sma, -largo_pin_exterior], 
                                                       [x_pos - desplazamiento_feed_x + radio_interior_sma, y_pos + radio_interior_sma, 0], priority=2)

                                inicio_puerto = [x_pos - desplazamiento_feed_x, y_pos, 0]
                                fin_puerto  = [x_pos - desplazamiento_feed_x, y_pos, espesor_sustrato]
                                
                                # Solo inyecta 1V al puerto activo, el resto 0V (Pasivos)
                                excitacion_voltaje = 1.0 if i == puerto_activo else 0.0
                                p = FDTD.AddLumpedPort(i+1, resistencia_alimentacion, inicio_puerto, fin_puerto, 'z', excitacion_voltaje, priority=5, edges2grid='xy')
                                puertos.append(p)

                            sustrato = CSX.AddMaterial('sustrato', epsilon=epsR_sustrato, kappa=kappa)
                            sustrato.AddBox([-ancho_sustrato/2, -largo_sustrato/2, 0], [ancho_sustrato/2, largo_sustrato/2, espesor_sustrato], priority=0)

                            plano_masa = CSX.AddMetal('plano_masa')
                            plano_masa.AddBox([-ancho_sustrato/2, -largo_sustrato/2, 0], [ancho_sustrato/2, largo_sustrato/2, 0], priority=10)

                            FDTD.AddEdges2Grid(dirs='xy', properties=parche, metal_edge_res=mesh_res/2)
                            mesh.SmoothMeshLines('all', mesh_res, 1.4)

                            nf2ff = FDTD.CreateNF2FFBox()
                            
                            archivo_csx = os.path.join(Sim_Path_Puerto, 'simp_patch_array.xml')
                            CSX.Write2XML(archivo_csx)

                            if puerto_activo == 0 and verificar_primera_geometria and primera_simulacion:
                                print("\n>>> ABRIENDO AppCSXCAD PARA VERIFICACIÓN VISUAL...")
                                os.system(f'AppCSXCAD "{archivo_csx}"')
                                respuesta = input("¿Continuar barrido automático? (s/n): ")
                                if respuesta.lower() != 's': import sys; sys.exit()
                                primera_simulacion = False

                            if not solo_post_proceso:
                                FDTD.Run(Sim_Path_Puerto, verbose=0, cleanup=True)

                            # 5. POST-PROCESO (FÍSICA DE HARDWARE)

                            if puerto_activo == 0:
                                f = np.linspace(max(1e9, f0-fc), f0+fc, 401)
                                puertos[0].CalcPort(Sim_Path_Puerto, f)
                                
                                s11_dB = 20.0 * np.log10(np.abs(puertos[0].uf_ref / (puertos[0].uf_inc + 1e-15)))
                                f_res_idx = np.argmin(s11_dB)
                                s11_minimo = s11_dB[f_res_idx]
                                f_res_calculada = f[f_res_idx] # FIJAMOS LA FRECUENCIA DE RESONANCIA PARA TODOS

                                # --- NUEVO: Cálculo de SWR ---
                                gamma = 10**(s11_minimo / 20.0)
                                swr_minimo = (1 + np.abs(gamma)) / (1 - np.abs(gamma))

                                # --- NUEVO: Cálculo de Cross Talk (Acoplamiento Mutuo) ---
                                acoplamientos = []
                                for p_idx in range(1, 5):
                                    puertos[p_idx].CalcPort(Sim_Path_Puerto, f)
                                    s_ix_dB = 20.0 * np.log10(np.abs(puertos[p_idx].uf_ref / (puertos[0].uf_inc + 1e-15)))
                                    acoplamientos.append(s_ix_dB[f_res_idx]) # Extrae acoplamiento en la f_res exacta
                                max_cross_talk = np.max(acoplamientos) # Peor caso de acoplamiento

                                plt.figure(figsize=(8,4))
                                plt.plot(f/1e9, s11_dB, 'k-')
                                plt.grid(True); plt.ylabel('S11 (dB)'); plt.xlabel('Frecuencia (GHz)')
                                plt.savefig(os.path.join(Sim_Path, '1_S11_Central.png')); plt.close()

                                # Far-Field del Puerto 0 (Hardware FoV)
                                nf2ff_res_0 = nf2ff.CalcNF2FF(Sim_Path_Puerto, f_res_calculada, theta_doa, [0., 90.], center=[0,0,0])
                                matriz_campos = nf2ff_res_0.E_norm[0]
                                
                                E_plane_norm = np.squeeze(20.0*np.log10(matriz_campos[:, 0] / (np.max(matriz_campos[:, 0]) + 1e-15)))
                                H_plane_norm = np.squeeze(20.0*np.log10(matriz_campos[:, 1] / (np.max(matriz_campos[:, 1]) + 1e-15)))
                                
                                hpbw_E = calcular_hpbw(theta_doa, E_plane_norm, caida_haz_db)
                                hpbw_H = calcular_hpbw(theta_doa, H_plane_norm, caida_haz_db)

                                plt.figure(figsize=(8,4))
                                plt.plot(theta_doa, E_plane_norm + nf2ff_res_0.Dmax[0], 'k-', label=f'E-plane FoV: {hpbw_E:.1f}°')
                                plt.plot(theta_doa, H_plane_norm + nf2ff_res_0.Dmax[0], 'r--', label=f'H-plane FoV: {hpbw_H:.1f}°')
                                plt.axhline(nf2ff_res_0.Dmax[0] - caida_haz_db, color='gray', linestyle=':')
                                plt.grid(True); plt.legend(); plt.ylabel('Directividad (dBi)'); plt.xlabel('Ángulo Theta (Grados)')
                                plt.savefig(os.path.join(Sim_Path, '2_FoV_Fisico.png'))
                                plt.close() # Cierra la figura para que no se imprima en consola

                            # Extracción estricta de Fase Compleja referenciada a [0,0,0]
                            nf2ff_calib = nf2ff.CalcNF2FF(Sim_Path_Puerto, f_res_calculada, theta_doa, [0.], center=[0,0,0])
                            manifold_complex[puerto_activo, :] = nf2ff_calib.E_theta[0][0, :] # Componente compleja Co-polar

                        # 6. POST-PROCESO DIGITAL (RESOLUCIÓN DOA)

                        print("    -> Calculando (MUSIC / Bartlett) ...")
                        
                        # NORMALIZACIÓN DEL MANIFOLD (Para que el SNR y el ruido sintético funcionen bien)
                        manifold_complex = manifold_complex / np.max(np.abs(manifold_complex))
                        
                        angulos_escenario = [-25.0, 25.0]
                        P_music, P_bartlett = simular_y_procesar_doa(manifold_complex, theta_doa, angulos_escenario, SNR_dB=15)
                        
                        res_bartlett = medir_ancho_pico_doa(theta_doa, P_bartlett, angulos_escenario[1])
                        res_music = medir_ancho_pico_doa(theta_doa, P_music, angulos_escenario[1])
                        rmse_music = calcular_rmse(theta_doa, P_music, angulos_escenario)

                        # Grafico del espectro de software
                        plt.figure(figsize=(10,5))
                        plt.plot(theta_doa, P_bartlett, 'b--', label=f'Bartlett (Res: {res_bartlett:.1f}°)')
                        plt.plot(theta_doa, P_music, 'r-', lw=2, label=f'MUSIC (Res: {res_music:.1f}°)')
                        plt.axvline(-10, color='k', linestyle=':', alpha=0.5)
                        plt.axvline(10, color='k', linestyle=':', alpha=0.5)
                        plt.title(f'Test de Resolución DOA (Cenit a $\pm10^\circ$)')
                        plt.xlabel('Ángulo Theta (Grados)'); plt.ylabel('Pseudo-Espectro (dB)')
                        plt.grid(True); plt.legend(); plt.ylim([-40, 5])
                        plt.savefig(os.path.join(Sim_Path, '3_Espectro_DOA.png')); plt.close()

                        # 7. REPORTE FINAL

                        res_B_str = f"{res_bartlett:.1f}" if res_bartlett != 999.0 else "FAIL"
                        res_M_str = f"{res_music:.1f}" if res_music != 999.0 else "FAIL"
                        rmse_M_str = f"{rmse_music:.2f}" if rmse_music != 999.0 else "FAIL"

                        config = {
                            'f0_GHz': float(f0/1e9), 'D_mm': float(D), 'L_mm': float(L_corte),
                            'Divisor': float(lambda_divisor), 'h_mm': float(espesor_sustrato),
                            'f_res_GHz': float(round(f_res_calculada/1e9, 3)), 'S11_dB': float(round(s11_minimo, 2)),
                            'SWR': float(round(swr_minimo, 2)), 'CT_dB': float(round(max_cross_talk, 2)),
                            'FoV_E': float(round(hpbw_E, 1)), 'Res_B': res_B_str, 'Res_M': res_M_str, 'RMSE': rmse_M_str
                        }
                        resultados_finales.append(config)
                        
                        with open(archivo_reporte, "a") as f:
                            f.write(f"{config['f0_GHz']:<9.3f} | {config['D_mm']:<7.1f} | {config['L_mm']:<7.1f} | {config['Divisor']:<5.1f} | {config['h_mm']:<7.1f} | {config['f_res_GHz']:<12.3f} | {config['S11_dB']:<9.2f} | {config['SWR']:<5.2f} | {config['CT_dB']:<8.2f} | {config['FoV_E']:<9.1f} | {config['Res_B']:<9} | {config['Res_M']:<9} | {config['RMSE']:<9}\n")

# 8. IMPRESIÓN FINAL EN CONSOLA
print("\n" + "="*138)
print(f"{'f0 (GHz)':<9} | {'D (mm)':<7} | {'L (mm)':<7} | {'Div':<5} | {'h (mm)':<7} | {'f_res (GHz)':<12} | {'S11 (dB)':<9} | {'SWR':<5} | {'CT (dB)':<8} | {'FoV_E(°)':<9} | {'Res_B(°)':<9} | {'Res_M(°)':<9} | {'RMSE(°)':<9}")
print("-" * 138)
for r in resultados_finales:
    print(f"{r['f0_GHz']:<9.3f} | {r['D_mm']:<7.1f} | {r['L_mm']:<7.1f} | {r['Divisor']:<5.1f} | {r['h_mm']:<7.1f} | {r['f_res_GHz']:<12.3f} | {r['S11_dB']:<9.2f} | {r['SWR']:<5.2f} | {r['CT_dB']:<8.2f} | {r['FoV_E']:<9.1f} | {r['Res_B']:<9} | {r['Res_M']:<9} | {r['RMSE']:<9}")
print("="*138 + "\n")
print(f"Todo guardado en: {carpeta_barrido}")